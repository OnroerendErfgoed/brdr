import logging
import math
import os.path

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
from brdr.constants import MULTI_SINGLE_ID_SEPARATOR
from brdr.constants import PERIMETER_ATTRIBUTE
from brdr.constants import SHAPE_INDEX_ATTRIBUTE
from brdr.constants import STABILITY
from brdr.constants import ZERO_STREAK
from brdr.enums import DiffMetric
from brdr.geometry_utils import buffer_neg
from brdr.geometry_utils import buffer_pos
from brdr.geometry_utils import from_crs
from brdr.geometry_utils import geometric_equality
from brdr.geometry_utils import get_bbox
from brdr.geometry_utils import get_geoms_from_geometry
from brdr.geometry_utils import get_partitions
from brdr.geometry_utils import get_shape_index
from brdr.geometry_utils import safe_intersection
from brdr.geometry_utils import safe_unary_union
from brdr.geometry_utils import to_crs
from brdr.geometry_utils import total_vertex_distance
from brdr.logger import LOGGER
from brdr.typings import ProcessResult

log = logging.getLogger(__name__)


def get_dict_geojsons_from_series_dict(
    series_dict: dict[str|int, dict[float, ProcessResult]],
    crs: CRS,
    id_field: str,
    series_prop_dict: dict[str|int, dict[float, str|int]] = None,
    geom_attributes=True,
):
    """
    Convert a series of process results to a GeoJSON feature collection.

    Args:
        series_dict (dict): Dictionary containing process results.
        crs (str): Coordinate reference system.
        id_field (str): Field name for the ID.
        series_prop_dict (dict, optional): Dictionary containing series properties.
        geom_attributes (bool, optional): Whether to include geometry attributes.

    Returns:
        dict: Dictionary of GeoJSON feature collections.
    """
    features_list_dict = {}

    for theme_id, results_dict in series_dict.items():

        prop_dict = dict(series_prop_dict or {}).get(theme_id, {})
        for relevant_distance, process_result in results_dict.items():
            properties = prop_dict.get(relevant_distance, {})
            properties[id_field] = theme_id
            for results_type, geom in process_result.items():
                if not isinstance(geom, BaseGeometry):
                    continue
                if results_type not in features_list_dict:
                    features_list_dict[results_type] = []

                feature = _feature_from_geom(
                    geom, theme_id, properties, geom_attributes
                )
                features_list_dict[results_type].append(feature)

    crs_geojson = {"type": "name", "properties": {"name": from_crs(crs)}}
    result = {
        result_type: FeatureCollection(features, crs=crs_geojson)
        for result_type, features in features_list_dict.items()
    }
    return result


def _feature_from_geom(
    geom: BaseGeometry,
    feature_id: str|int,
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


def geojson_from_dict(dictionary, crs: CRS, id_field, prop_dict=None, geom_attributes=True):
    """
    Get a GeoJSON (FeatureCollection) from a dictionary of IDs (keys) and geometries (values).

    Args:
        dictionary (dict): Dictionary of geometries.
        crs (str): Coordinate reference system.
        id_field (str): Field name for the ID.
        prop_dict (dict, optional): Dictionary of properties.
        geom_attributes (bool, optional): Whether to include geometry attributes.

    Returns:
        FeatureCollection: The GeoJSON FeatureCollection.
    """
    features = []
    for key, geom in dictionary.items():
        properties = dict(prop_dict or {}).get(key, {})
        properties[id_field] = key
        features.append(_feature_from_geom(geom, key, properties, geom_attributes))
    crs_geojson = {"type": "name", "properties": {"name": from_crs(crs)}}
    geojson = FeatureCollection(features, crs=crs_geojson)
    return geojson


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


def multi_to_singles(dict_geoms):
    """
    Convert a dictionary of Shapely (multi-)geometries to a dictionary containing only single geometries.

    Args:
        dict_geoms (dict): Dictionary of geometries.

    Returns:
        tuple: A tuple containing:
            - dict: Dictionary of single geometries.
            - dict: Dictionary mapping new keys to original keys.
    """
    resulting_dict_geoms = {}
    dict_multi_as_single = {}
    for key, geom in dict_geoms.items():
        if geom is None or geom.is_empty:
            continue

        geometries = list(get_geoms_from_geometry(geom))
        if len(geometries) == 1:
            resulting_dict_geoms[key] = geometries[0]
            continue

        i = 0
        for g in geometries:
            new_key = str(key) + MULTI_SINGLE_ID_SEPARATOR + str(i)
            dict_multi_as_single[new_key] = key
            resulting_dict_geoms[new_key] = g
            i = i + 1
    return resulting_dict_geoms, dict_multi_as_single


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


def diffs_from_dict_processresult(
    dict_processresult: dict[float, ProcessResult],
    geom_thematic: BaseGeometry,
    reference_union: BaseGeometry,
    diff_metric: DiffMetric = DiffMetric.CHANGES_AREA,
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
        diff = diff_from_processresult(
            processresult, geom_thematic, reference_union, diff_metric
        )

        diffs[rel_dist] = diff
    return diffs


def diff_from_processresult(processresult, geom_thematic, reference_union, diff_metric):
    if geom_thematic.geom_type in (
        "LineString",
        "MultiLineString",
    ):
        diff_metric = DiffMetric.CHANGES_LENGTH
        # diff_metric = DiffMetric.REFERENCE_USAGE
        # diff_metric = DiffMetric.TOTAL_DISTANCE
    elif geom_thematic.geom_type in (
        "Point",
        "MultiPoint",
    ):
        diff_metric = DiffMetric.TOTAL_DISTANCE
        diff_metric = DiffMetric.REFERENCE_USAGE
        diff_metric = DiffMetric.TOTAL_DISTANCE
    result = processresult.get("result")
    result_diff = processresult.get("result_diff")
    diff = 0
    original = geom_thematic
    if result_diff is None or result_diff.is_empty or result is None or result.is_empty:
        diff = 0
    elif diff_metric == DiffMetric.TOTAL_AREA:
        diff = result.area - original.area
    elif diff_metric == DiffMetric.TOTAL_PERCENTAGE:
        diff = result.area - original.area
        diff = diff * 100 / result.area
    elif diff_metric == DiffMetric.CHANGES_AREA:
        diff = result_diff.area
    elif diff_metric == DiffMetric.CHANGES_PERCENTAGE:
        diff = result_diff.area
        diff = diff * 100 / result.area
    elif diff_metric == DiffMetric.TOTAL_LENGTH:
        diff = result.length - original.length
    elif diff_metric == DiffMetric.CHANGES_LENGTH:
        diff = result_diff.length
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
    diff = round(diff, 1)
    return diff


def fetch_all_ogc_features(base_url, params=None, headers=None, max_pages=math.inf):
    """
    Fetches all features from an OGC Feature API, regardless of the pagination mechanism.

    Supports:
    - Cursor-based pagination via 'next' links
    - Offset-based pagination via 'startIndex'
    - No pagination (all features in a single page)

    :param base_url: URL of the /items endpoint of the OGC API
    :param headers: Optional headers (e.g., Accept: application/json)
    :return: A list of all features
    """

    if params is None:
        params = {}
    if headers is None:
        headers = {"Accept": "application/json"}

    all_features = []
    url = base_url
    limit = params.get("limit", DOWNLOAD_LIMIT)
    start_index = 0
    local_params = params.copy()
    use_start_index = False
    page = 0

    while url and page <= max_pages:
        LOGGER.debug("page:" + str(page) + "url: " + str(url))
        # Add startIndex when using startIndex
        if use_start_index:
            local_params["startIndex"] = start_index

        response = requests.get(url, params=local_params, headers=headers)
        response.raise_for_status()
        data = response.json()

        # Add features
        all_features.extend(data.get("features", []))

        # Search next-link
        next_link = next(
            (link["href"] for link in data.get("links", []) if link["rel"] == "next"),
            None,
        )

        if next_link:
            # Cursor-based pagination
            url = next_link
            local_params = None  # next-url has all its parameters included
        elif (
            "numberMatched" in data
            and "numberReturned" in data
            and data["numberMatched"] == data["numberReturned"]
        ):
            break
        elif "numberReturned" in data and "features" in data:
            # limit not used and all features returned
            if data["numberReturned"] > limit:
                break
            # Possibly offset-based pagination
            if not use_start_index:
                use_start_index = True
                start_index = params.get("startIndex", 0)
                local_params = params.copy()
                local_params["count"] = limit
            start_index += limit
            if data["numberReturned"] <= limit:
                break
        else:
            # No further pagination
            break

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


def geojson_to_dicts(collection, id_property):
    """
    Converts a GeoJSON collection into dictionaries of geometries and properties.

    Parameters:
    collection (dict): The GeoJSON collection to convert.
    id_property (str): The property name to use as the key for the dictionaries.

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
        key = f["properties"][
            id_property
        ]  # TODO to check if this has to be converted to string?
        geom = shape(f["geometry"])
        data_dict[key] = make_valid(geom)
        data_dict_properties[key] = f["properties"]
    return data_dict, data_dict_properties


def get_collection_by_partition(
    url,
    params,
    geometry,
    partition=1000,
    crs=DEFAULT_CRS,
):
    """
    Retrieves a collection of geographic data by partitioning the input geometry.

    Parameters:
    url (str): The base URL for the data source.
    geometry (object): The geometric area to partition and retrieve data for. If None, retrieves data for the entire area.
    partition (int, optional): The number of partitions to divide the geometry into. Default is 1000. If less than 1, no partitioning is done.
    limit (int, optional): The maximum number of items to retrieve. Default is DOWNLOAD_LIMIT.
    crs (str, optional): The coordinate reference system to use. Default is DEFAULT_CRS.

    Returns:
    dict: A collection of geographic data, potentially partitioned by the input geometry.
    """
    crs = to_crs(crs)
    collection = {}
    if geometry is None or geometry.is_empty:
        collection = get_collection(url=url, params=params)
    elif partition < 1:
        params["bbox"] = get_bbox(geometry)
        params["bbox-crs"] = from_crs(crs)
        collection = get_collection(url=url, params=params)
    else:
        geoms = get_partitions(geometry, partition)
        for g in geoms:
            params["bbox"] = get_bbox(g)
            params["bbox-crs"] = from_crs(crs)
            coll = get_collection(url=url, params=params)
            if collection == {}:
                collection = dict(coll)
            elif "features" in collection and "features" in coll:
                collection["features"].extend(coll["features"])
    return collection


def merge_process_results(
    result_dict: dict[str|int, dict[float, ProcessResult]], dict_multi_as_single: dict
) -> dict[str|int, dict[float, ProcessResult]]:
    """
     Merges processresults in a dictionary from multiple themeIDs into a single themeID.

    Args: result_dict (dict): A dictionary where keys are theme IDs and values are
        process results

    Returns: dict: A new dictionary with merged geometries and remarks (processresults), where keys are global
        theme IDs and values are merged geometries and remarks.

    """
    grouped_results: dict[str|int, dict[float, ProcessResult]] = {}

    for id_theme, dict_results in result_dict.items():
        if id_theme in dict_multi_as_single.keys():
            id_theme_global = dict_multi_as_single[id_theme]
        else:
            id_theme_global = id_theme
        if id_theme_global not in grouped_results:
            grouped_results[id_theme_global] = dict_results
        else:
            for rel_dist, process_result in dict_results.items():
                for key in process_result:
                    value = process_result[key]  # noqa
                    if isinstance(value, str) and value != "":
                        existing_remark: str = grouped_results[id_theme_global][
                            rel_dist
                        ][
                            key
                        ]  # noqa
                        grouped_results[id_theme_global][rel_dist][key] = (
                            existing_remark + " | " + str(value)
                        )
                        continue
                    elif isinstance(value, BaseGeometry):
                        geom = value
                        if geom.is_empty or geom is None:
                            continue
                        existing: BaseGeometry = grouped_results[id_theme_global][
                            rel_dist
                        ][
                            key
                        ]  # noqa
                        grouped_results[id_theme_global][rel_dist][key] = (
                            safe_unary_union([existing, geom])
                        )  # noqa
    return grouped_results


def is_brdr_formula(brdr_formula):
    """
    returns true if the value has the correct structure of a base_formula, otherwise False
    :param brdr_formula:
    :return:
    """
    if brdr_formula is None or not isinstance(brdr_formula, dict):
        return False
    if brdr_formula.keys() >= {
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


def unary_union_result_dict(result_dict):
    # make a unary union for each key value in the result dict
    for key in ProcessResult.__annotations__:
        geom = result_dict.get(key, GeometryCollection())  # noqa
        if isinstance(geom, BaseGeometry) and not geom.is_empty:
            geom = safe_unary_union(geom)
        result_dict[key] = geom  # noqa
    return result_dict


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


def equal_geom_in_array(geom, geom_array, correction_distance, mitre_limit):
    """
    Check if a predicted geometry is equal to other predicted geometries in a list.
    Equality is defined as there is the symmetrical difference is smaller than the CORRECTION DISTANCE
    Returns True if one of the elements is equal, otherwise False
    """
    for g in geom_array:
        if geometric_equality(geom, g, correction_distance, mitre_limit):
            return True
    return False
