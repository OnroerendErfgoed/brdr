import hashlib
import json
import logging
import math
import os.path
import uuid
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import requests
from geojson import Feature
from geojson import FeatureCollection
from geojson import dump
from pyproj import CRS
from shapely import GeometryCollection
from shapely import make_valid
from shapely import node
from shapely import polygonize
from shapely.geometry import shape
from shapely.geometry.base import BaseGeometry

from brdr.constants import AREA_ATTRIBUTE
from brdr.constants import DEFAULT_CRS
from brdr.constants import DOWNLOAD_LIMIT
from brdr.constants import PERIMETER_ATTRIBUTE
from brdr.constants import SHAPE_INDEX_ATTRIBUTE
from brdr.constants import STABILITY
from brdr.constants import ZERO_STREAK
from brdr.enums import DiffMetric
from brdr.geometry_utils import buffer_neg
from brdr.geometry_utils import buffer_pos
from brdr.geometry_utils import from_crs
from brdr.geometry_utils import get_bbox
from brdr.geometry_utils import get_partitions
from brdr.geometry_utils import get_shape_index
from brdr.geometry_utils import safe_intersection
from brdr.geometry_utils import safe_unary_union
from brdr.geometry_utils import to_crs
from brdr.geometry_utils import total_vertex_distance
from brdr.logger import LOGGER
from brdr.typings import ProcessResult

log = logging.getLogger(__name__)


def get_geojsons_from_process_results(
    process_results: dict[str | int, dict[float, ProcessResult]],
    crs: CRS,
    id_field: str,
    series_prop_dict: dict[str | int, dict[float, str | int]] = None,
    geom_attributes=True,
):
    """
    Convert a series of process results to a GeoJSON feature collection.

    Args:
        process_results (dict): Dictionary containing process results.
        crs (str): Coordinate reference system.
        id_field (str): Field name for the ID.
        series_prop_dict (dict, optional): Dictionary containing series properties.
        geom_attributes (bool, optional): Whether to include geometry attributes.

    Returns:
        dict: Dictionary of GeoJSON feature collections.
    """
    features_list_dict = {}

    for theme_id, results_dict in process_results.items():

        prop_dict = dict(series_prop_dict or {}).get(theme_id, {})
        for relevant_distance, process_result in results_dict.items():
            properties = prop_dict.get(relevant_distance, {})
            properties[id_field] = theme_id
            properties.update(process_result["properties"])
            for results_type, value in process_result.items():
                if not isinstance(value, BaseGeometry):
                    continue
                if results_type not in features_list_dict:
                    features_list_dict[results_type] = []

                feature = _feature_from_geom(
                    value, theme_id, properties, geom_attributes
                )
                features_list_dict[results_type].append(feature)

    crs_geojson = {"type": "name", "properties": {"name": from_crs(crs)}}
    geojsons = {
        result_type: FeatureCollection(features, crs=crs_geojson)
        for result_type, features in features_list_dict.items()
    }
    return geojsons


def deep_merge(dict_a, dict_b):
    for key, value in dict_b.items():
        if key in dict_a and isinstance(dict_a[key], dict) and isinstance(value, dict):
            deep_merge(dict_a[key], value)
        else:
            dict_a[key] = value
    return dict_a


def _feature_from_geom(
    geom: BaseGeometry,
    feature_id: str | int,
    properties: dict = None,
    geom_attributes=True,
) -> Feature:
    """
    Convert a geometry to a GeoJSON feature.

    Args:
        geom (BaseGeometry): The geometry to convert.
        feature_id: a (unique) id for the feature
        properties (dict, optional): The properties to include in the feature.
        geom_attributes (bool, optional): Whether to include geometry attributes.

    Returns:
        Feature: The GeoJSON feature.
    """
    properties = dict(properties or {})
    if geom_attributes:
        area = geom.area
        perimeter = geom.length
        properties[AREA_ATTRIBUTE] = area
        properties[PERIMETER_ATTRIBUTE] = perimeter
        properties[SHAPE_INDEX_ATTRIBUTE] = get_shape_index(area, perimeter)
    return Feature(geometry=geom, id=feature_id, properties=properties)


def write_geojson(path_to_file, geojson):
    """
    Write a GeoJSON object to a file.

    Args:
        path_to_file (str): Path to the output file.
        geojson (FeatureCollection): The GeoJSON object to write.
    """
    parent = os.path.dirname(path_to_file)
    os.makedirs(parent, exist_ok=True)
    with open(path_to_file, "w") as f:
        dump(geojson, f, default=str)


def flatten_iter(lst):
    while any(isinstance(i, list) for i in lst):
        lst = [
            item
            for sublist in lst
            for item in (sublist if isinstance(sublist, list) else [sublist])
        ]
    return lst


def polygonize_reference_data(dict_ref):
    """
    Create a new dictionary with non-overlapping polygons based on a reference data dictionary.

    Args:
        dict_ref (dict): Dictionary of reference geometries.

    Returns:
        dict: Dictionary of non-overlapping polygons.
    """
    arr_ref = []
    for key in dict_ref:
        arr_ref.append((make_valid(dict_ref[key])))
    arr_ref = node(GeometryCollection(arr_ref))
    polygonized = polygonize(arr_ref.geoms)
    dict_ref = {}
    i = 0
    for g in polygonized.geoms:
        if str(g.geom_type) in ["Polygon", "MultiPolygon"]:
            i = i + 1
            key = str(i)  # unique keys will be added (reference_id is lost)
            dict_ref[key] = make_valid(g)
        else:
            logging.warning("geom excluded: " + str(g))
    return dict_ref


def coverage_ratio(values, min_val=0, bin_count=10):
    max_val = max(values)
    if len(values) == 0 or max_val <= 0:
        return 0.0
    bin_size = round((max_val - min_val) / bin_count, 2)
    bins = np.arange(min_val, max_val + bin_size, bin_size)

    # Tel hoeveel bins minstens één waarde bevatten
    bin_counts, _ = np.histogram(values, bins)
    filled_bins = np.count_nonzero(bin_counts)

    total_bins = len(bins) - 1  # histogram geeft n-1 bins
    return filled_bins / total_bins


def determine_stability(x, y):
    """
    Determine the stability and indicators ('zero_streaks') based on the derivative of the difference measurement (y) of a xy-plot

    Args:
        x (numpy.ndarray): The x values of the plot: relevant distances
        y (numpy.ndarray): The y values of the plot: diffence measurement

    Returns:
        tuple: A tuple containing:
            - list: List of breakpoints (extremes).
            - list: List of zero_streaks.
    """
    derivative = _numerical_derivative(x, y)
    max_zero_streak_score = x[-1] - x[0]
    start_streak = None
    streak = 0
    write_zero_streak = False
    dict_stability = dict()
    for i in range(1, len(x)):
        dict_stability[x[i - 1]] = dict()
        dict_stability[x[i - 1]][STABILITY] = False
        dict_stability[x[i - 1]][ZERO_STREAK] = None
        if round(derivative[i], 2) == 0:  # noqa
            dict_stability[x[i - 1]][STABILITY] = True
            streak = streak + 1
            if start_streak is None:
                start_streak = x[i - 1]
        elif streak != 0:
            write_zero_streak = True
        if start_streak is not None and (write_zero_streak or len(x) - 1 == i):
            end_streak = x[i - 1]
            if derivative[i] == 0:
                end_streak = x[i]
            score_streak = round(
                (end_streak - start_streak) * 100 / max_zero_streak_score, 2
            )
            center_streak = start_streak + (end_streak - start_streak) / 2

            dict_stability[start_streak][ZERO_STREAK] = (
                start_streak,
                end_streak,
                center_streak,
                score_streak,
            )

            streak = 0
            start_streak = None
            write_zero_streak = False
            logging.debug("end_streak")

    return dict_stability


def _numerical_derivative(x, y):
    """
    Calculate the numerical derivative of a graph.

    Parameters:
      x (numpy.ndarray): The x values of the graph.
      y (numpy.ndarray): The y values of the graph.

    Returns:
      numpy.ndarray: The derivative of y with respect to x.
    """

    dx = np.diff(x)
    dy = np.diff(y)
    derivative = dy / dx
    derivative = np.insert(derivative, 0, 0)
    # derivative = np.append(derivative, 0)

    return derivative


def get_geometry_difference_metrics_from_processresults(
    dict_processresult: dict[float, ProcessResult],
    geom_thematic: BaseGeometry,
    reference_union: BaseGeometry,
    diff_metric: DiffMetric = DiffMetric.SYMMETRICAL_AREA_CHANGE,
):
    """
    Calculates a dictionary containing difference metrics for thematic geometry based on a distance series.

    Parameters:
    dict_processresult (dict): A dictionary where keys are thematic IDs and values are dictionaries mapping relevant distances to ProcessResult objects.
    geom_thematic (BaseGeometry): A dictionary where keys are thematic IDs and values are BaseGeometry objects representing the original geometries.
    reference_union (BaseGeometry): unioned reference geometries
    diff_metric (DiffMetric, optional): The metric to use for calculating differences. Default is DiffMetric.CHANGES_AREA.

    Returns:
    dict: A dictionary where keys are relevant distances to calculated difference metrics.
    """

    diffs = {}

    for rel_dist in dict_processresult:
        processresult = dict_processresult.get(rel_dist, {})
        if geom_thematic.geom_type in (
            "LineString",
            "MultiLineString",
        ):
            diff_metric = DiffMetric.LENGTH_CHANGE
        elif geom_thematic.geom_type in (
            "Point",
            "MultiPoint",
        ):
            diff_metric = DiffMetric.TOTAL_DISTANCE
        diff = get_geometry_difference_metrics_from_processresult(
            processresult, geom_thematic, reference_union, diff_metric
        )

        diffs[rel_dist] = diff
    return diffs


def get_geometry_difference_metrics_from_processresult(
    processresult, geom_thematic, reference_union, diff_metric
):
    result = processresult.get("result")
    result_diff = processresult.get("result_diff")
    result_diff_min = processresult.get("result_diff_min")
    diff = 0
    original = geom_thematic
    if result_diff is None or result_diff.is_empty or result is None or result.is_empty:
        diff = 0
    elif diff_metric == DiffMetric.SYMMETRICAL_AREA_CHANGE:
        diff = result_diff.area
    elif diff_metric == DiffMetric.SYMMETRICAL_AREA_PERCENTAGE_CHANGE:
        diff = result_diff.area
        if original.area != 0:
            diff = diff * 100 / original.area
    elif diff_metric == DiffMetric.AREA_CHANGE:
        diff = result.area - original.area
    elif diff_metric == DiffMetric.AREA_PERCENTAGE_CHANGE:
        diff = result.area - original.area
        if original.area != 0:
            diff = diff * 100 / original.area
    elif diff_metric == DiffMetric.LENGTH_CHANGE:
        diff = result.length - original.length
    elif diff_metric == DiffMetric.LENGTH_PERCENTAGE_CHANGE:
        diff = result.length - original.length
        if original.length != 0:
            diff = diff * 100 / original.length
    elif diff_metric == DiffMetric.LENGTH_REMOVED:
        if not result_diff_min is None or result_diff_min.is_empty:
            diff = result_diff_min.length
    elif diff_metric == DiffMetric.REFERENCE_USAGE:
        if not reference_union is None and not reference_union.is_empty:
            reference_union_buffer = buffer_pos(reference_union, 0.01)
            result_buffer = buffer_pos(result, 0.01)
            reference_usage_geom = safe_intersection(
                result_buffer, reference_union_buffer
            )
        else:
            reference_usage_geom = None
        if reference_usage_geom is not None and not reference_usage_geom.is_empty:
            diff = reference_usage_geom.area
        else:
            diff = 0
    elif diff_metric == DiffMetric.TOTAL_DISTANCE:
        diff = total_vertex_distance(original, result, bidirectional=False)

    # round, so the detected changes are within 10cm, 10cm² or 0.1%
    diff = round(abs(diff), 1)
    return diff


def fetch_all_ogc_features(
    base_url, params=None, headers=None, max_pages=math.inf, max_workers=None
):
    """
    Fetches all features from an OGC Feature API using parallel requests where possible.

    If the API provides 'numberMatched', the function calculates all offsets and
    fetches pages in parallel. If not, it falls back to sequential cursor-based
    pagination via 'next' links.

    Parameters
    ----------
    base_url : str
        URL of the /items endpoint of the OGC API.
    params : dict, optional
        Query parameters for the API request. Default is None.
    headers : dict, optional
        Optional headers (e.g., Accept: application/json).
    max_pages : int or math.inf, optional
        Maximum number of pages to fetch. Default is math.inf.
    max_workers : int, optional
        Maximum number of concurrent threads for parallel fetching. Default is None.

    Returns
    -------
    list
        A list containing all features retrieved from the API.
    """
    if params is None:
        params = {}
    if headers is None:
        headers = {"Accept": "application/json"}

    limit = params.get("limit", DOWNLOAD_LIMIT)

    # 1. Initial request to determine total count and pagination type
    response = requests.get(base_url, params=params, headers=headers)
    response.raise_for_status()
    data = response.json()

    all_features = list(data.get("features", []))
    total_matched = data.get("numberMatched")
    number_returned = data.get("numberReturned", len(data.get("features", [])))

    # Check for 'next' link for cursor-based fallback
    next_link = next(
        (link["href"] for link in data.get("links", []) if link["rel"] == "next"),
        None,
    )

    # 2. Parallel Path: If we know the total and it's more than one page
    if total_matched and total_matched > number_returned and max_pages > 0:
        LOGGER.info(
            f"Total features: {total_matched}. Switching to parallel offset-based retrieval."
        )

        # Calculate offsets
        offsets = range(number_returned, total_matched, limit)
        # Limit offsets by max_pages
        offsets = list(offsets)[: int(max_pages)]

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_offset = {}
            for offset in offsets:
                local_params = params.copy()
                local_params["startIndex"] = offset
                local_params["limit"] = limit
                future_to_offset[
                    executor.submit(
                        requests.get, base_url, params=local_params, headers=headers
                    )
                ] = offset

            for future in as_completed(future_to_offset):
                try:
                    resp = future.result()
                    resp.raise_for_status()
                    page_data = resp.json()
                    all_features.extend(page_data.get("features", []))
                except Exception as exc:
                    LOGGER.error(
                        f"Offset {future_to_offset[future]} generated an exception: {exc}"
                    )

        return all_features

    # 3. Sequential Path: Fallback for Cursor-based pagination
    elif next_link and max_pages > 0:
        LOGGER.info(
            "Parallel retrieval not possible (no total count). Falling back to sequential cursor pagination."
        )
        url = next_link
        page = 1
        while url and page <= max_pages:
            LOGGER.debug(f"Fetching page {page}: {url}")
            # For next_links, parameters are usually already in the URL
            response = requests.get(url, headers=headers)
            response.raise_for_status()
            data = response.json()

            all_features.extend(data.get("features", []))

            url = next(
                (
                    link["href"]
                    for link in data.get("links", [])
                    if link["rel"] == "next"
                ),
                None,
            )
            page += 1

    return all_features


def make_feature_collection(features):
    """
    Maakt een GeoJSON FeatureCollection van een lijst met features.

    :param features: lijst van GeoJSON features (dicts met 'type': 'Feature')
    :return: dict met 'type': 'FeatureCollection' en 'features': [...]
    """
    return {"type": "FeatureCollection", "features": features}


def get_collection(url, params=None):
    """
    Fetches a collection of features from a (paginated) API endpoint (OGC Feature API or WFS).

    This function retrieves a collection of features from a URL that supports
    pagination using a `startIndex` parameter. It iteratively retrieves features in
    chunks of the specified `limit` until no more features are available.

    Args:
        url (str): The base URL of the API endpoint.
        params: parameters for the request

    Returns:
        dict: A dictionary representing the complete GeoJSON feature collection.
            This might be truncated if the total number of features exceeds the
            limitations of the API or server.
    """
    if params is None:
        params = {}
    features = fetch_all_ogc_features(base_url=url, params=params)
    return make_feature_collection(features)


def geojson_geometry_to_shapely(geojson_geometry):
    """
    Converts a geojson geometry into a shapely geometry
    :param geojson_geometry:  geojson geometry
    :return: shapely geometry
    """
    return shape(geojson_geometry)


def geojson_to_dicts(collection, id_property=None):
    """
    Converts a GeoJSON collection into dictionaries of geometries and properties.

    Parameters:
    collection (dict): The GeoJSON collection to convert.
    id_property (str): The property name to use as the key for the dictionaries. if None, the feature-identifier 'id' is used when available.

    Returns:
    tuple: Two dictionaries:
        - data_dict (dict): A dictionary where keys are the id_property values and values are the geometries.
        - data_dict_properties (dict): A dictionary where keys are the id_property values and values are the properties.
    """
    data_dict = {}
    data_dict_properties = {}
    if collection is None or "features" not in collection:
        return data_dict, data_dict_properties
    for f in collection["features"]:
        if not id_property:
            if "id" in f:
                key = f["id"]
            else:
                raise KeyError(
                    "Feature-identifier 'id' not found in GeoJson FeatureCollection. Please provide a Geojson with feature-identifiers or define id_property (name of attribute with unique ids) from the properties"
                )
        else:
            if "properties" in f and id_property in f["properties"]:
                key = f["properties"][id_property]
            else:
                raise KeyError(
                    f"Identifier '{id_property}' not found in properties of GeoJson FeatureCollection"
                )
        geom = shape(f["geometry"])
        data_dict[key] = make_valid(geom)
        data_dict_properties[key] = f["properties"]
    return data_dict, data_dict_properties


def get_collection_by_partition(
    url, params, geometry, partition=1000, crs=DEFAULT_CRS, max_workers=None
):
    """
    Retrieves a collection of geographic data by partitioning the input geometry.

    This function parallelizes data retrieval by splitting the geometry into
    smaller parts and fetching each partition concurrently using threads.

    Parameters
    ----------
    url : str
        The base URL for the OGC API Features data source.
    params : dict
        The query parameters for the API request.
    geometry : object
        The geometric area (Shapely object) to partition and retrieve data for.
        If None or empty, retrieves data for the entire area.
    partition : int, optional
        The number of partitions to divide the geometry into. Default is 1000.
        If less than 1, no partitioning is performed.
    crs : str, optional
        The coordinate reference system to use. Default is DEFAULT_CRS.
    max_workers : int, optional
        The maximum number of concurrent threads to use for parallel requests.
        Default is None.

    Returns
    -------
    dict
        A GeoJSON-like collection of geographic data containing merged features
        from all partitions.
    """

    crs = to_crs(crs)
    collection = {}

    if geometry is None or geometry.is_empty:
        return get_collection(url=url, params=params)

    if partition < 1:
        local_params = params.copy()
        local_params["bbox"] = get_bbox(geometry)
        local_params["bbox-crs"] = from_crs(crs)
        return get_collection(url=url, params=local_params)

    # Prepare partitions
    geoms = get_partitions(geometry, partition)
    features_list = []

    # Use ThreadPoolExecutor for parallel I/O-bound API requests
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_geom = {}

        for g in geoms:
            # Create a shallow copy of params to avoid race conditions
            thread_params = params.copy()
            thread_params["bbox"] = get_bbox(g)
            thread_params["bbox-crs"] = from_crs(crs)

            # Submit the request to the thread pool
            future = executor.submit(get_collection, url=url, params=thread_params)
            future_to_geom[future] = g

        # Collect results as they complete
        for future in as_completed(future_to_geom):
            try:
                coll = future.result()
                if coll and "features" in coll:
                    # If this is the first valid response, use it as the template for the collection
                    if not collection:
                        collection = dict(coll)
                        collection["features"] = []

                    features_list.extend(coll["features"])
            except Exception as exc:
                # Retrieve the specific geometry that caused the failure for debugging
                geom_info = future_to_geom[future]
                print(f"Partition request generated an exception: {exc}")
                raise exc

    # Re-attach all collected features to the template
    if collection:
        unique_features = deduplicate_features(features_list)
        collection["features"] = unique_features
    # TODO; only return fatures that intersect with the 'geometry'?

    return collection


def deduplicate_features(features):
    """
    Removes duplicate features using ID or a content-based hash as fallback.

    Parameters
    ----------
    features : list
        List of GeoJSON-like feature dictionaries.

    Returns
    -------
    list
        Deduplicated list of features.
    """
    seen_hashes = set()
    unique_features = []

    for feat in features:
        # 1. Try to get the official ID
        feat_id = feat.get("id")

        # 2. If no ID, create a hash of the content
        if feat_id is None:
            # We use a stable JSON string to represent the feature content
            content = json.dumps(
                {"g": feat.get("geometry"), "p": feat.get("properties")}, sort_keys=True
            )
            feat_id = hashlib.md5(content.encode()).hexdigest()

        if feat_id not in seen_hashes:
            unique_features.append(feat)
            seen_hashes.add(feat_id)

    return unique_features


def build_reverse_index_wkb(d: dict):
    return {g.wkb: k for k, g in d.items()}


def is_brdr_observation(brdr_observation):
    """
    returns true if the value has the correct structure of a base_observation, otherwise False
    :param brdr_observation:
    :return:
    """
    if brdr_observation is None or not isinstance(brdr_observation, dict):
        return False
    if brdr_observation.keys() >= {
        "alignment_date",
        "brdr_version",
        "reference_source",
        "full",
        "area",
        "reference_features",
        "reference_od",
    }:
        return True
    return False


def union_process_result(process_result: ProcessResult) -> ProcessResult:
    """
    Perform a unary union on all geometries within a ProcessResult object.

    This function iterates through all annotated geometric fields in the
    ProcessResult, attempts to merge overlapping or adjacent geometries
    into a single unified structure using a safe unary union operation,
    and updates the result in place.

    Parameters
    ----------
    process_result : ProcessResult
        A dictionary-like or dataclass object containing geometric data
        (typically resulting from a polygon processing operation).
        The keys correspond to specific geometric categories or layers.

    Returns
    -------
    ProcessResult
        The updated ProcessResult object where each valid, non-empty
        geometry has been replaced by its unioned representation.

    Notes
    -----
    - Empty geometries or those that are not instances of `BaseGeometry`
      are ignored.
    - The function utilizes `safe_unary_union` to handle potential
      topological inconsistencies or invalid geometries that might
      cause a standard union to fail.
    - Fields are accessed based on the type annotations of the
      `ProcessResult` class.

    See Also
    --------
    safe_unary_union : The underlying function used for geometric merging.
    ProcessorConfig : Settings that define how these results are generated.
    """
    for key in ProcessResult.__annotations__:
        if key in ["properties", "metadata", "observation"]:
            process_result[key] = process_result.get(key, {})  # noqa
            continue
        value = process_result.get(key, GeometryCollection())  # noqa
        if isinstance(value, BaseGeometry):
            value = safe_unary_union(value)
        process_result[key] = value  # noqa

    return process_result


def get_relevant_polygons_from_geom(
    geometry: BaseGeometry, buffer_distance: float, mitre_limit
):
    """
    Get only the relevant parts (polygon) from a geometry.
    Points, Lines and Polygons smaller than relevant distance are excluded from the
    result
    """
    if not geometry or geometry.is_empty:
        # If the input geometry is empty or None, do nothing.
        return geometry
    else:
        geometry = safe_unary_union(geometry)
        # Create a GeometryCollection from the input geometry.
        geometry_collection = GeometryCollection(geometry)
        array = []
        for g in geometry_collection.geoms:
            # Ensure each sub-geometry is valid.
            g = make_valid(g)
            if str(g.geom_type) in ["Polygon", "MultiPolygon"]:
                relevant_geom = buffer_neg(
                    g,
                    buffer_distance,
                    mitre_limit=mitre_limit,
                )
                if relevant_geom is not None and not relevant_geom.is_empty:
                    array.append(g)
    return safe_unary_union(array)


def urn_from_geom(geom: BaseGeometry):
    return uuid.UUID(hex=hashlib.sha256(geom.wkb).hexdigest()[::2]).urn


# def equal_geom_in_array(geom, geom_array, correction_distance, mitre_limit):
#     """
#     Check if a predicted geometry is equal to other predicted geometries in a list.
#     Equality is defined as there is the symmetrical difference is smaller than the CORRECTION DISTANCE
#     Returns True if one of the elements is equal, otherwise False
#     """
#     for g in geom_array:
#         if geometric_equality(geom, g, correction_distance, mitre_limit):
#             return True
#     return False


# def create_full_interpolated_dataset(
#     discrete_list: List[float],  # The complete list of all x-coordinates
#     cached_results: Dict[
#         float, Any
#     ],  # The dictionary with the calculated (cached) values
# ) -> Dict[float, Any]:
#     """
#     Populates the complete discrete list by interpolating non-calculated values
#     with the nearest calculated boundary value (Zero-Order Hold).
#
#     Returns: A new dictionary with results for ALL x-coordinates in the list.
#     """
#
#     full_results = {}
#     x_coords = np.array(discrete_list)
#
#     # We start with the value of the first calculated point
#     current_fill_value = None
#
#     # Find the very first calculated point to determine the starting value
#     for x in x_coords:
#         if x in cached_results:
#             # Make a deepcopy of the very first value found.
#             current_fill_value = copy.deepcopy(cached_results[x])
#             break
#
#     if current_fill_value is None:
#         # This shouldn't happen if the list has boundaries, but for safety
#         return full_results
#
#     # 1. Iterate through the entire discrete list and perform the interpolation
#     for x in x_coords:
#         if x in cached_results:
#             # 2. Point is calculated (this is an interval boundary in the stable region).
#             # Replace the fill_value with a DEEPCOPY of the new cached value.
#             current_fill_value = copy.deepcopy(cached_results[x])
#             # Assign a DEEPCOPY of the new fill_value to full_results.
#             full_results[x] = copy.deepcopy(current_fill_value)
#         else:
#             # 3. Point is NOT calculated (lies within a stable region).
#             # We populate it with the value of the last calculated (left boundary) point.
#             # Assign a DEEPCOPY of the fill_value to full_results.
#             full_results[x] = copy.deepcopy(current_fill_value)
#
#     return full_results


# def recursive_stepwise_interval_check(
#     f_value: Callable[
#         [float], Any
#     ],  # Function that returns the NUMERIC value/Object (V)
#     f_condition: Callable[
#         [Any, Any], bool
#     ],  # Function that performs the BOOLEAN assessment on V_start and V_end
#     discrete_list: List[float],  # The FIXED, discrete list of x-coordinates
#     initial_sample_size: int = None,  # Number of evenly spaced values to initially cache
# ) -> Tuple[List[float], Dict[float, Any]]:
#     """
#     Recursively searches a discrete domain (x-coordinates) to find intervals where a
#     stability condition (f_condition) based on calculated values (f_value) is violated.
#
#     It uses a Zero-Order Hold (ZOH) approach for condition checking and a cache-aware
#     stability check to prevent premature recursion stop in unstable regions.
#
#     The check is performed by evaluating the f_condition between interval endpoints
#     and recursively subdividing the interval if the condition is not met.
#
#     Args:
#         f_value: A function that calculates and returns the complex/numerical
#                  result (V) for a given x-coordinate.
#         f_condition: A function that takes V_start and V_end and returns True if
#                      the stability/convergence condition is met, False otherwise.
#         discrete_list: The complete, sorted list of x-coordinates to be evaluated.
#         initial_sample_size: Optional. If set, this number of points are pre-calculated
#                              and cached to potentially speed up the initial checks.
#
#     Returns:
#         A tuple containing:
#         1. A sorted list of x-coordinates that were identified as boundaries of
#            the non-stable (non-matching) regions.
#         2. The complete dictionary cache of all calculated f_value results.
#     """
#
#     x_coords = np.array(discrete_list)
#     cache: Dict[float, Any] = {}
#     non_matching_x = set()
#
#     def _get_f_value(x: float) -> Any:
#         """Helper to retrieve f_value result from cache or calculate it."""
#         if x not in cache:
#             # print(f"  > F_VALUE(x) calculated for x={round(x, 4)}")
#             # We use deepcopy here to ensure the cached object is independent if f_value
#             # returns a mutable object that might be modified later outside this function.
#             # However, typically f_value is expected to return a stable result object.
#             # If the result itself is used as input for further modification *within*
#             # the surrounding logic, deepcopy may be added, but based on the original
#             # context, we assume a stable return value for now.
#             cache[x] = f_value(x)
#         return cache[x]
#
#     # --- Initial Sample (Pre-caching) ---
#     if (
#         initial_sample_size is not None
#         and initial_sample_size > 0
#         and initial_sample_size < len(x_coords)
#     ):
#         step = max(1, len(x_coords) // initial_sample_size)
#         initial_indices = np.arange(0, len(x_coords), step)[:initial_sample_size]
#         for x in x_coords[initial_indices]:
#             _get_f_value(x)
#         # print(f"✅ Initialization: {len(initial_indices)} points cached.")
#     # -----------------------------------
#
#     def _recursive_check(start_index: int, end_index: int):
#         """
#         The core recursive function to check the condition in the interval [start_index, end_index].
#         """
#
#         # 0. Base Case: Interval contains no intermediate values (adjacent points)
#         if end_index - start_index <= 1:
#             x_start, x_end = x_coords[start_index], x_coords[end_index]
#             num_start, num_end = _get_f_value(x_start), _get_f_value(x_end)
#
#             if not f_condition(num_start, num_end):
#                 non_matching_x.add(x_start)
#                 non_matching_x.add(x_end)
#             return
#
#         x_start, x_end = x_coords[start_index], x_coords[end_index]
#         num_start, num_end = _get_f_value(x_start), _get_f_value(x_end)
#
#         interval_is_satisfied = f_condition(num_start, num_end)
#
#         needs_forced_recursion = False
#
#         if interval_is_satisfied:
#             # **NEW LOGIC:** Check internal stability using already cached points
#
#             # Find all x-values in the cache that fall STRICTLY within this interval
#             internal_cached_x = sorted([x for x in cache.keys() if x_start < x < x_end])
#
#             # Iterate through the sub-intervals [current_x, next_x]
#             current_x_val = x_start
#             current_num_val = num_start
#
#             # Combine the internal cached points with the endpoint to cover all sub-intervals
#             all_check_points = internal_cached_x + [x_end]
#
#             for next_x_val in all_check_points:
#                 # We know next_x_val is either in the cache or is x_end (which is already calculated)
#                 next_num_val = cache[next_x_val]
#
#                 # Evaluate the sub-interval [current_x_val, next_x_val]
#                 if not f_condition(current_num_val, next_num_val):
#                     # Internal instability found! Override the lazy stop.
#                     needs_forced_recursion = True
#                     break
#
#                 current_x_val = next_x_val
#                 current_num_val = next_num_val
#
#         # 3. Decision Point:
#
#         if not interval_is_satisfied or needs_forced_recursion:
#             # Condition NOT satisfied OR Internal cache indicates instability: Continue recursion.
#
#             # Add the extreme values (they are part of the 'non-matching' region)
#             non_matching_x.add(x_start)
#             non_matching_x.add(x_end)
#
#             # Perform recursion on the halves (Stepwise refinement)
#             mid_index = (start_index + end_index) // 2
#
#             _recursive_check(start_index, mid_index)
#             _recursive_check(mid_index, end_index)
#
#         else:
#             # Condition IS satisfied and NO internal problems detected: Stop recursion.
#             return
#
#     # Start the recursion over the full discrete list
#     _recursive_check(0, len(x_coords) - 1)
#
#     return sorted(list(non_matching_x)), cache
