import logging
import os.path

import numpy as np
import requests
from geojson import Feature
from geojson import FeatureCollection
from geojson import dump
from shapely import GeometryCollection
from shapely import Polygon
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
from brdr.geometry_utils import buffer_neg_pos
from brdr.geometry_utils import buffer_pos
from brdr.geometry_utils import geometric_equality
from brdr.geometry_utils import get_bbox
from brdr.geometry_utils import get_geoms_from_geometry
from brdr.geometry_utils import get_partitions
from brdr.geometry_utils import get_shape_index
from brdr.geometry_utils import safe_difference
from brdr.geometry_utils import safe_intersection
from brdr.geometry_utils import safe_unary_union
from brdr.geometry_utils import snap_geometry_to_reference
from brdr.geometry_utils import total_vertex_distance
from brdr.typings import ProcessResult

log = logging.getLogger(__name__)


def get_dict_geojsons_from_series_dict(
    series_dict: dict[any, dict[float, ProcessResult]],
    crs: str,
    id_field: str,
    series_prop_dict: dict[any, dict[float, any]] = None,
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

    crs_geojson = {"type": "name", "properties": {"name": crs}}
    result = {
        result_type: FeatureCollection(features, crs=crs_geojson)
        for result_type, features in features_list_dict.items()
    }
    return result


def _feature_from_geom(
    geom: BaseGeometry,
    feature_id: any,
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


def geojson_from_dict(dictionary, crs, id_field, prop_dict=None, geom_attributes=True):
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
    crs_geojson = {"type": "name", "properties": {"name": crs}}
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


def get_collection(ref_url, limit):
    """
    Fetches a collection of features from a paginated API endpoint.

    This function retrieves a collection of features from a URL that supports
    pagination using a `startIndex` parameter. It iteratively retrieves features in
    chunks of the specified `limit` until no more features are available.

    Args:
        ref_url (str): The base URL of the API endpoint.
        limit (int): The maximum number of features to retrieve per request.

    Returns:
        dict: A dictionary representing the complete GeoJSON feature collection.
            This might be truncated if the total number of features exceeds the
            limitations of the API or server.

    Logs:
        - Debug logs the URL being used for each request during pagination.
    """
    start_index = 0
    collection = {}
    while True:
        url = ref_url + "&startIndex=" + str(start_index)
        logging.debug("called url: " + url)
        json = requests.get(url).json()
        feature_collection = dict(json)
        if (
            "features" not in feature_collection
            or len(feature_collection["features"]) == 0
        ):
            break
        start_index = start_index + limit
        if collection == {}:
            collection = feature_collection
        else:
            collection["features"].extend(feature_collection["features"])  # noqa
        if len(feature_collection["features"]) < limit:
            break
    return collection


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
    url, geometry, partition=1000, limit=DOWNLOAD_LIMIT, crs=DEFAULT_CRS
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
    collection = {}
    if geometry is None or geometry.is_empty:
        collection = get_collection(
            _add_bbox_to_url(url=url, crs=crs, bbox=None), limit
        )
    elif partition < 1:
        collection = get_collection(
            _add_bbox_to_url(url=url, crs=crs, bbox=get_bbox(geometry)), limit
        )
    else:
        geoms = get_partitions(geometry, partition)
        for g in geoms:
            coll = get_collection(
                _add_bbox_to_url(url=url, crs=crs, bbox=get_bbox(g)), limit
            )
            if collection == {}:
                collection = dict(coll)
            elif "features" in collection and "features" in coll:
                collection["features"].extend(coll["features"])
    return collection


def _add_bbox_to_url(url, crs=DEFAULT_CRS, bbox=None):
    """
    Adds a bounding box (bbox) parameter to the URL for geographic data requests.

    Parameters:
    url (str): The base URL for the data source.
    crs (str, optional): The coordinate reference system to use. Default is DEFAULT_CRS.
    bbox (str, optional): The bounding box coordinates to add to the URL. If None, no bbox is added.

    Returns:
    str: The updated URL with the bbox parameter included, if provided.
    """
    # Load the Base reference data
    if bbox is not None:
        url = url + "&bbox-crs=" + crs + "&bbox=" + bbox
    return url


def merge_process_results(
    result_dict: dict[any, dict[float, ProcessResult]], dict_multi_as_single: dict
) -> dict[any, dict[float, ProcessResult]]:
    """
     Merges processresults in a dictionary from multiple themeIDs into a single themeID.

    Args: result_dict (dict): A dictionary where keys are theme IDs and values are
        process results

    Returns: dict: A new dictionary with merged geometries and remarks (processresults), where keys are global
        theme IDs and values are merged geometries and remarks.

    """
    grouped_results: dict[any, dict[float, ProcessResult]] = {}

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


def _unary_union_result_dict(result_dict):
    # make a unary union for each key value in the result dict
    for key in ProcessResult.__annotations__:
        geom = result_dict.get(key, GeometryCollection())  # noqa
        if isinstance(geom, BaseGeometry) and not geom.is_empty:
            geom = safe_unary_union(geom)
        result_dict[key] = geom  # noqa
    return result_dict


def _calculate_geom_by_intersection_and_reference(
    geom_intersection: BaseGeometry,
    geom_reference: BaseGeometry,
    input_geometry_inner: BaseGeometry,
    is_open_domain,
    buffer_distance,
    threshold_overlap_percentage,
    threshold_exclusion_percentage,
    threshold_exclusion_area,
    threshold_inclusion_percentage,
    mitre_limit,
    partial_snapping,
    partial_snap_strategy,
    partial_snap_max_segment_length,
):
    """
    Calculates the geometry based on intersection and reference geometries.

    Args:
        geom_intersection (BaseGeometry): The intersection geometry.
        geom_reference (BaseGeometry): The reference geometry.
        is_open_domain (bool): A flag indicating whether it's a public domain
            (area not covered with reference polygon).
        threshold_exclusion_percentage (int): The threshold exclusion percentage.
        threshold_exclusion_area (int): The threshold exclusion area.
        buffer_distance (float): The buffer distance.
        threshold_overlap_percentage (int): The threshold overlap percentage.

    Returns:
        tuple: A tuple containing the resulting geometries:

        *   geom: BaseGeometry or None: The resulting geometry or None if conditions
            are not met.
        *   geom_relevant_intersection: BaseGeometry or None: The relevant
            intersection.
        *   geom_relevant_difference: BaseGeometry or None: The relevant difference.

    Notes:
        -   If the reference geometry area is 0, the overlap is set to 100%.
        -   If the overlap is less than relevant_OVERLAP_PERCENTAGE or the
            intersection area is less than relevant_OVERLAP_AREA, None is returned.
        -   Otherwise, the relevant intersection and difference geometries are
            calculated.
        -   If both relevant intersection and difference are non-empty, the final
            geometry is obtained by applying safe intersection and buffering.
        -   If only relevant intersection is non-empty, the result is the reference
            geometry.
        -   If only relevant difference is non-empty, the result is None.
    """
    od_overlap = 111  # define a specific value for defining overlap of OD
    if geom_reference.area == 0:
        overlap = od_overlap  # Open Domain

    else:
        overlap = geom_intersection.area * 100 / geom_reference.area

    if (
        overlap < threshold_exclusion_percentage
        or geom_intersection.area < threshold_exclusion_area
    ):
        return Polygon(), Polygon(), Polygon()

    if overlap >= threshold_inclusion_percentage and not overlap == od_overlap:
        return geom_reference, geom_reference, Polygon()

    geom_difference = safe_difference(geom_reference, geom_intersection)
    geom_relevant_intersection = buffer_neg(
        geom_intersection,
        buffer_distance,
        mitre_limit=mitre_limit,
    )
    geom_intersection_inner = safe_intersection(
        geom_intersection, input_geometry_inner
    )  # this part is the intersection-part that always has to be kept, because it is inside the inner_geometry

    geom_relevant_difference = buffer_neg(
        geom_difference,
        buffer_distance,
        mitre_limit=mitre_limit,
    )

    if (
        not geom_intersection_inner.is_empty
        and geom_relevant_intersection.is_empty
        and not geom_relevant_difference.is_empty
    ):
        geom_x = safe_intersection(
            geom_intersection, buffer_pos(geom_intersection_inner, 2 * buffer_distance)
        )

        geom_x = snap_geometry_to_reference(
            geom_x,
            geom_reference,
            max_segment_length=partial_snap_max_segment_length,
            snap_strategy=partial_snap_strategy,
            tolerance=2 * buffer_distance,
        )

        geom = geom_x
    elif (
        not geom_relevant_intersection.is_empty
        and not geom_relevant_difference.is_empty
    ):
        # relevant intersection and relevant difference

        geom_x = safe_difference(
            geom_reference,
            safe_intersection(
                geom_difference,
                buffer_neg_pos(
                    geom_difference,
                    buffer_distance,
                    mitre_limit=mitre_limit,
                ),
            ),
        )
        geom_x = buffer_neg_pos(geom_x, buffer_distance, mitre_limit=mitre_limit)

        geom_intersection_buffered = buffer_pos(geom_intersection, 2 * buffer_distance)
        geom_difference_2 = safe_difference(geom_reference, geom_intersection_buffered)
        geom_difference_2_buffered = buffer_pos(geom_difference_2, 2 * buffer_distance)

        geom_x = safe_difference(geom_x, geom_difference_2_buffered)

        geom_x = safe_intersection(geom_x, geom_reference)

        if partial_snapping:
            geom_x = snap_geometry_to_reference(
                geom_x,
                geom_reference,
                max_segment_length=partial_snap_max_segment_length,
                snap_strategy=partial_snap_strategy,
                tolerance=2 * buffer_distance,
            )
        geom = safe_unary_union(
            [geom_x, geom_relevant_intersection, geom_intersection_inner]
        )

        # when calculating for OD, we create a 'virtual parcel'. When calculating this
        # virtual parcel, it is buffered to take outer boundaries into account.
        # This results in a side effect that there are extra non-logical parts included
        # in the result. The function below tries to exclude these non-logical parts.
        # see eo_id 206363 with relevant distance=0.2m and SNAP_ALL_SIDE
        if is_open_domain:
            geom = _get_relevant_polygons_from_geom(geom, buffer_distance, mitre_limit)
    elif not geom_relevant_intersection.is_empty and geom_relevant_difference.is_empty:
        geom = geom_reference
    elif geom_relevant_intersection.is_empty and not geom_relevant_difference.is_empty:
        geom = geom_relevant_intersection  # (=empty geometry)
    else:
        # No relevant intersection and no relevant difference
        if is_open_domain:
            # geom = geom_relevant_intersection  # (=empty geometry)
            geom = snap_geometry_to_reference(
                geom_intersection,
                geom_reference,
                snap_strategy=partial_snap_strategy,
                tolerance=2 * buffer_distance,
                max_segment_length=partial_snap_max_segment_length,
            )
        elif not geom_intersection_inner.is_empty:
            geom_intersection_buffered = buffer_pos(
                geom_intersection, 2 * buffer_distance
            )
            geom_difference_2 = safe_difference(
                geom_reference, geom_intersection_buffered
            )
            geom_difference_2_buffered = buffer_pos(
                geom_difference_2, 2 * buffer_distance
            )
            geom = safe_difference(geom_reference, geom_difference_2_buffered)
        elif threshold_overlap_percentage < 0:
            # if we take a value of -1, the original border will be used
            geom = geom_intersection
        elif overlap > threshold_overlap_percentage:
            geom = geom_reference
        else:
            geom = geom_relevant_intersection  # (=empty geometry)
    return geom, geom_relevant_intersection, geom_relevant_difference


def _get_relevant_polygons_from_geom(
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


def _equal_geom_in_array(geom, geom_array, correction_distance, mitre_limit):
    """
    Check if a predicted geometry is equal to other predicted geometries in a list.
    Equality is defined as there is the symmetrical difference is smaller than the CORRECTION DISTANCE
    Returns True if one of the elements is equal, otherwise False
    """
    for g in geom_array:
        if geometric_equality(geom, g, correction_distance, mitre_limit):
            return True
    return False
