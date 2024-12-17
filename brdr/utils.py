import logging
import os.path

import numpy as np
import requests
from geojson import Feature, FeatureCollection, dump
from shapely import GeometryCollection, make_valid, node, polygonize, unary_union
from shapely.geometry import shape
from shapely.geometry.base import BaseGeometry

from brdr.constants import (
    MULTI_SINGLE_ID_SEPARATOR,
    DEFAULT_CRS,
    DOWNLOAD_LIMIT,
    RELEVANT_DISTANCE_FIELD_NAME,
    NR_CALCULATION_FIELD_NAME,
    REMARK_FIELD_NAME,
    PERIMETER_ATTRIBUTE,
    SHAPE_INDEX_ATTRIBUTE,
    AREA_ATTRIBUTE,
)
from brdr.enums import DiffMetric
from brdr.geometry_utils import get_partitions, get_bbox, get_shape_index
from brdr.typings import ProcessResult


def get_series_geojson_dict(
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
        nr_calculations = len(results_dict)
        prop_dict = dict(series_prop_dict or {}).get(theme_id, {})
        for relative_distance, process_result in results_dict.items():
            properties = prop_dict.get(relative_distance, {})
            properties[id_field] = theme_id
            properties[NR_CALCULATION_FIELD_NAME] = nr_calculations
            properties[RELEVANT_DISTANCE_FIELD_NAME] = relative_distance
            if "remark" in process_result:
                properties[REMARK_FIELD_NAME] = process_result["remark"]

            for results_type, geom in process_result.items():
                if not isinstance(geom, BaseGeometry):
                    continue
                if results_type not in features_list_dict:
                    features_list_dict[results_type] = []

                feature = _feature_from_geom(geom, properties, geom_attributes)
                features_list_dict[results_type].append(feature)

    crs_geojson = {"type": "name", "properties": {"name": crs}}
    return {
        result_type: FeatureCollection(features, crs=crs_geojson)
        for result_type, features in features_list_dict.items()
    }


def _feature_from_geom(
    geom: BaseGeometry,
    properties: dict = None,
    geom_attributes=True,
) -> Feature:
    """
    Convert a geometry to a GeoJSON feature.

    Args:
        geom (BaseGeometry): The geometry to convert.
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
    return Feature(geometry=geom, properties=properties)


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
        features.append(_feature_from_geom(geom, properties, geom_attributes))
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


def multipolygons_to_singles(dict_geoms):
    """
    Convert a dictionary of Shapely geometries to a dictionary containing only single polygons.

    Args:
        dict_geoms (dict): Dictionary of geometries.

    Returns:
        tuple: A tuple containing:
            - dict: Dictionary of single polygons.
            - dict: Dictionary mapping new keys to original keys.
    """
    resulting_dict_geoms = {}
    dict_multi_as_single = {}
    for key, geom in dict_geoms.items():
        if geom is None or geom.is_empty:
            continue
        elif str(geom.geom_type) == "Polygon":
            resulting_dict_geoms[key] = geom
        elif str(geom.geom_type) == "MultiPolygon":
            polygons = list(geom.geoms)
            if len(polygons) == 1:
                resulting_dict_geoms[key] = polygons[0]
                continue
            i = 0
            for p in polygons:
                new_key = str(key) + MULTI_SINGLE_ID_SEPARATOR + str(i)
                dict_multi_as_single[new_key] = key
                resulting_dict_geoms[new_key] = p
                i = i + 1
        else:
            logging.debug("geom excluded: " + str(geom) + " for key: " + str(key))
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


def get_breakpoints_zerostreak(x, y):
    """
    Determine the extremes and zero_streaks of a graph based on the derivative.

    Args:
        x (numpy.ndarray): The x values of the graph.
        y (numpy.ndarray): The y values of the graph.

    Returns:
        tuple: A tuple containing:
            - list: List of breakpoints (extremes).
            - list: List of zero_streaks.
    """
    derivative = _numerical_derivative(x, y)
    # plt.plot(x, y, label="y")
    # plt.plot( x, derivative, label="derivative")
    # plt.legend()
    # plt.show()
    extremes = []
    zero_streaks = []
    start_streak = None
    streak = 0
    write_zero_streak = False
    last_extreme = 0
    for i in range(1, len(x)):
        if round(derivative[i], 2) == 0:  # noqa
            streak = streak + 1
            if start_streak is None:
                start_streak = x[i - 1]
        elif streak != 0:
            write_zero_streak = True
        if start_streak is not None and (write_zero_streak or len(x) - 1 == i):
            end_streak = x[i - 1]
            if derivative[i] == 0:
                end_streak = x[i]
            center_streak = start_streak + (end_streak - start_streak) / 2
            zero_streaks.append(
                (start_streak, end_streak, center_streak, streak, last_extreme)
            )
            streak = 0
            start_streak = None
            write_zero_streak = False
            logging.debug("end_streak")
        if (
            i < len(x) - 1
            and round(derivative[i], 2) > 0  # noqa
            and derivative[i - 1] <= derivative[i]
            and derivative[i] >= derivative[i + 1]
        ):
            last_extreme = derivative[i]
            extremes.append((x[i], derivative[i], "maximum"))
        if (
            i < len(x) - 1
            and round(derivative[i], 2) < 0  # noqa
            and derivative[i - 1] >= derivative[i]
            and derivative[i] <= derivative[i + 1]
        ):
            last_extreme = derivative[i]
            extremes.append((x[i], derivative[i], "minimum"))
    for extremum in extremes:
        logging.debug(
            f"breakpoints: relevant_distance:"
            f"{extremum[0]:.2f}, extreme:{extremum[1]:.2f} ({extremum[2]})"
        )
    for st in zero_streaks:
        logging.debug(
            f"zero_streaks: [{st[0]:.2f} - {st[1]:.2f}] - center:{st[2]:.2f}"
            f" - counter:{st[3]:.2f} - min/max-extreme:{st[4]:.2f} "
        )
    # plt.plot(series, afgeleide, label='afgeleide-' + str(key))
    return extremes, zero_streaks


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


def diffs_from_dict_series(
    dict_series: dict[any, dict[float, ProcessResult]],
    dict_thematic: dict[any, BaseGeometry],
    diff_metric: DiffMetric = DiffMetric.CHANGES_AREA,
):
    """
    Calculates a dictionary containing difference metrics for thematic elements based on a distance series.

    Parameters:
    dict_series (dict): A dictionary where keys are thematic IDs and values are dictionaries mapping relative distances to ProcessResult objects.
    dict_thematic (dict): A dictionary where keys are thematic IDs and values are BaseGeometry objects representing the original geometries.
    diff_metric (DiffMetric, optional): The metric to use for calculating differences. Default is DiffMetric.CHANGES_AREA.

    Returns:
    dict: A dictionary where keys are thematic IDs and values are dictionaries mapping relative distances to calculated difference metrics.
    """
    diffs = {}
    # all the relevant distances used to calculate the series
    for thematic_id, results_dict in dict_series.items():
        diffs[thematic_id] = {}

        for rel_dist in results_dict:
            result = results_dict.get(rel_dist, {}).get("result")
            result_diff = results_dict.get(rel_dist, {}).get("result_diff")

            diff = 0
            if (
                result_diff is None
                or result_diff.is_empty
                or result is None
                or result.is_empty
            ):
                diff = 0
            elif diff_metric == DiffMetric.TOTAL_AREA:
                diff = result.area - dict_thematic[thematic_id].area
            elif diff_metric == DiffMetric.TOTAL_PERCENTAGE:
                diff = result.area - dict_thematic[thematic_id].area
                diff = diff * 100 / result.area
            elif diff_metric == DiffMetric.CHANGES_AREA:
                # equals the symmetrical difference, so equal to
                # result_diff_plus.area + result_diff_min.area
                diff = result_diff.area
                # diff = result_diff_plus.area + result_diff_min.area
            elif diff_metric == DiffMetric.CHANGES_PERCENTAGE:
                diff = result_diff.area
                diff = diff * 100 / result.area

            # round, so the detected changes are within 10cm² or 0.1%
            diff = round(diff, 1)
            diffs[thematic_id][rel_dist] = diff

    return diffs


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
        key = str(f["properties"][id_property])
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
                        grouped_results[id_theme_global][rel_dist][key] = unary_union(
                            [existing, geom]
                        )  # noqa
    return grouped_results
