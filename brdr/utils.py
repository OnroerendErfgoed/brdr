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
from shapely import unary_union
from shapely.geometry import shape
from shapely.geometry.base import BaseGeometry

from brdr.constants import MULTI_SINGLE_ID_SEPARATOR
from brdr.enums import DiffMetric
from brdr.typings import ProcessResult


# def geojson_tuple_from_tuple(
#     my_tuple, crs, name_id, prop_dict=None, geom_attributes=True
# ):
#     """
#     get a geojson-tuple (6 geojsons) for a tuple of results (results, result_diff, ...)
#     """
#     feature_collections = []
#     for count, tup in enumerate(my_tuple):
#         feature_collections.append(
#             geojson_from_dict(
#                 tup, crs, name_id, prop_dict=prop_dict, geom_attributes=geom_attributes
#             )
#         )
#     return tuple(feature_collections)


def get_series_geojson_dict(
    series_dict: dict[float, dict[str, ProcessResult]],
    crs: str,
    id_field: str,
    series_prop_dict: dict[float, dict[str, any]] = None,
    geom_attributes=True,
):
    """
    Convert a series of process results to a GeoJSON feature collections.
    """
    features_list_dict = {}

    for relative_distance, results_dict in series_dict.items():
        prop_dict = (series_prop_dict or {}).get(relative_distance, {})
        for theme_id, process_result in results_dict.items():
            properties = prop_dict.get(theme_id, {})
            properties[id_field] = theme_id
            properties["relevant_distance"] = relative_distance

            for results_type, geom in process_result.items():
                if results_type not in features_list_dict:
                    features_list_dict[results_type] = []

                feature = feature_from_geom(geom, properties, geom_attributes)
                features_list_dict[results_type].append(feature)

    crs_geojson = {"type": "name", "properties": {"name": crs}}
    return {
        result_type: FeatureCollection(features, crs=crs_geojson)
        for result_type, features in features_list_dict.items()
    }


def feature_from_geom(
    geom: BaseGeometry,
    properties: dict = None,
    geom_attributes=True,
) -> Feature:
    """
    Convert a geometry to a GeoJSON feature.

    Args:
        geom (BaseGeometry): The geometry to convert.
        properties (dict): The properties to include in the feature.
        geom_attributes (bool): Whether to include geometry attributes (default True).

    Returns:
        Feature: The GeoJSON feature.
    """

    properties = properties or {}
    if geom_attributes:
        area = geom.area
        perimeter = geom.length
        properties["area"] = area
        properties["perimeter"] = perimeter
        properties["shape_index"] = perimeter / area if area != 0 else -1
    return Feature(geometry=geom, properties=properties)


def geojson_from_dict(dictionary, crs, id_field, prop_dict=None, geom_attributes=True):
    """
    get a geojson (featurecollection) from a dictionary of ids(keys) and geometries (values)
    """
    features = []
    for key, geom in dictionary.items():
        properties = (prop_dict or {}).get(key, {})
        properties[id_field] = key
        features.append(feature_from_geom(geom, properties, geom_attributes))
    crs_geojson = {"type": "name", "properties": {"name": crs}}
    geojson = FeatureCollection(features, crs=crs_geojson)
    return geojson


def write_geojson(path_to_file, geojson):
    parent = os.path.dirname(path_to_file)
    os.makedirs(parent, exist_ok=True)
    with open(path_to_file, "w") as f:
        dump(geojson, f)


def multipolygons_to_singles(dict_geoms):
    """
    Converts a dictionary of shapely-geometries to a dictionary containing only single polygons.

    This function iterates through a dictionary where values are Shapely-geometries and performs the following:

    * **Polygons:** Preserves the key and geometry from the original dictionary.
    * **MultiPolygons with one polygon:** Preserves the key and extracts the single polygon.
    * **MultiPolygons with multiple polygons:**
        * Creates new keys for each polygon by appending a suffix (_index) to the original key.
        * Assigns the individual polygons from the MultiPolygon to the newly created keys.

    Args:
        dict_geoms (dict): A dictionary where keys are identifiers and values are GeoJSON geometries.

    Returns:
        dict: A new dictionary containing only single polygons (as Polygon geometries).
              Keys are created based on the logic described above.

    Notes:
        * Geometries that are not Polygons or MultiPolygons are excluded with a warning message printed.
    """
    resulting_dict_geoms = {}
    for key in dict_geoms:
        geom = dict_geoms[key]
        if str(geom.geom_type) == "Polygon":
            resulting_dict_geoms[key] = geom
        elif str(geom.geom_type) == "MultiPolygon":
            polygons = list(geom.geoms)
            if len(polygons) == 1:
                resulting_dict_geoms[key] = polygons[0]
                continue
            i = 0
            for p in polygons:
                new_key = str(key) + MULTI_SINGLE_ID_SEPARATOR + str(i)
                resulting_dict_geoms[new_key] = p
                i = i + 1
        else:
            print("geom excluded: " + str(geom))
    return resulting_dict_geoms


def polygonize_reference_data(dict_ref):
    """
    Creates a new dictionary with non-overlapping polygons based on a reference data dictionary.

    This function is designed to handle situations where the original reference data dictionary might contain:

    * Overlapping polygons: It creates new, non-overlapping polygons by combining all reference borders.
    * Multiple overlapping references: This function is useful when combining references like parcels and buildings that might overlap.

    **Important:** The original reference IDs are lost in the process of creating new non-overlapping polygons. New unique keys are assigned instead.

    Args:
        dict_ref (dict): A dictionary where keys are identifiers and values are Shapely geometries (assumed to be Polygons or MultiPolygons).

    Returns:
        dict: A new dictionary containing non-overlapping polygons derived from the original reference data.
              Keys are unique strings (reference IDs are lost).

    Notes:
        * Geometries that are not Polygons or MultiPolygons are excluded with a warning message printed.
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


def get_oe_dict_by_ids(objectids, oetype="aanduidingsobjecten"):
    """
    Fetches thematic data for a list of objectIDs from the Inventaris Onroerend Erfgoed API.

    This function retrieves information about designated heritage objects (erfgoedobjecten or aanduidingsobjecten)
    from the Flemish Agency for Heritage (Inventaris Onroerend Erfgoed) based on a list of their IDs.

    Args:
        objectids (list): A list of objectIDs of 'erfgoedobjecten' or 'aanduidingsobjecten'.
        oetype (string): A string: 'aanduidingsobjecten' (default) or 'erfgoedobjecten'

    Returns:
        dict: A dictionary where keys are objectIDs (as strings) and values are GeoJSON geometry objects.
              If an erfgoedobject/aanduidingsobject is not found, a corresponding warning message will be logged
              but it won't be included in the returned dictionary.

    Raises:
        requests.exceptions.RequestException: If there is an error fetching data from the API.
    """
    dict_thematic = {}
    base_url = "https://inventaris.onroerenderfgoed.be/" + oetype + "/"
    headers = {"Accept": "application/json"}
    for a in objectids:
        url = base_url + str(a)
        response = requests.get(url, headers=headers).json()
        if "id" in response.keys():
            key = str(response["id"])
            geom = shape(response["locatie"]["contour"])
            dict_thematic[key] = geom
        else:
            logging.warning("object id " + str(a) + " not available in " + oetype)
    return dict_thematic


def get_oe_geojson_by_bbox(bbox, limit=1000):
    """
    Fetches GeoJSON data for designated heritage objects (aanduidingsobjecten) within a bounding box.

    This function retrieves information about aanduidingsobjecten from the Flemish
    Mercator public WFS service using a bounding box (bbox) as a filter. The bbox should
    be provided in the format "xmin,ymin,xmax,ymax" (EPSG:31370 projection).

    Args:
        bbox (str): A comma-separated string representing the bounding box in EPSG:31370
                   projection (e.g., "100000,500000,200000,600000").
        limit (int, optional): The maximum number of features to retrieve per request.
                              Defaults to 1000.

    Returns:
        dict: A dictionary containing the retrieved GeoJSON feature collection. This
              collection might be truncated if the total number of features exceeds
              the specified limit.
    """
    theme_url = (
        "https://www.mercator.vlaanderen.be/raadpleegdienstenmercatorpubliek/wfs?"
        "SERVICE=WFS&REQUEST=GetFeature&VERSION=2.0.0&"
        "TYPENAMES=ps:ps_aandobj&"
        f"COUNT={str(limit)}&"
        "SRSNAME=urn:ogc:def:crs:EPSG::31370&"
        f"BBOX={bbox}&outputFormat=application/json"
    )
    start_index = 0
    collection = {}
    while True:
        url = theme_url + "&startIndex=" + str(start_index)
        feature_collection = requests.get(url).json()
        if (
            "features" not in feature_collection
            or len(feature_collection["features"]) == 0
        ):
            break
        start_index = start_index + limit
        collection = collection | feature_collection
    return collection


def get_breakpoints_zerostreak(x, y):
    """
    Determine the extremes and zero_streaks of a graph based on the derivative, and return:
    * the breakpoints: extremes (breakpoints) of graph where 'change' occurs
    * the zero_streaks: ranges where the derivative is zero, ranges of relevant_distance where 'no-change' occurs

    Parameters:
    x (numpy.ndarray): The x values of the graph.
    derivative (numpy.ndarray): The y values of the graph.

    Returns:
    extremes: A list of tuples for breakpoints:
        *   relevant distance where extreme occurs
        *   extreme value
        *   minimum or maximum
    zero_streaks: A list of tuples for zero_streaks:
        * relevant distance where zero_streak starts
        * relevant distance where zero_streak ends
        * center of start- and end- zero_streak
        * counter of #relevant_distances where zero-streak holds on
        * extreme value for zero_streak
    """
    derivative = numerical_derivative(x, y)
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
        if round(derivative[i], 2) == 0:
            streak = streak + 1
            if start_streak == None:
                start_streak = x[i - 1]
        elif streak != 0:
            write_zero_streak = True
        if start_streak != None and (write_zero_streak or len(x) - 1 == i):
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
            and round(derivative[i], 2) > 0
            and derivative[i - 1] <= derivative[i]
            and derivative[i] >= derivative[i + 1]
        ):
            last_extreme = derivative[i]
            extremes.append((x[i], derivative[i], "maximum"))
        if (
            i < len(x) - 1
            and round(derivative[i], 2) < 0
            and derivative[i - 1] >= derivative[i]
            and derivative[i] <= derivative[i + 1]
        ):
            last_extreme = derivative[i]
            extremes.append((x[i], derivative[i], "minimum"))
    for extremum in extremes:
        logging.info(
            f"breakpoints: relevant_distance:{extremum[0]:.2f}, extreme:{extremum[1]:.2f} ({extremum[2]})"
        )
    for st in zero_streaks:
        logging.info(
            f"zero_streaks: [{st[0]:.2f} - {st[1]:.2f}] - center:{st[2]:.2f} - counter:{st[3]:.2f} - min/max-extreme:{st[4]:.2f} "
        )
    # plt.plot(series, afgeleide, label='afgeleide-' + str(key))
    return extremes, zero_streaks


def numerical_derivative(x, y):
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


def _filter_dict_by_key(dictionary, filter_key):
    """
    Filters a dictionary to only include keys matching a specific value.

    This function creates a new dictionary containing entries from the original dictionary
    where the key matches the provided `filter_key`.

    Args:
        dictionary (dict): The dictionary to filter.
        filter_key (str): The key value to filter by.

    Returns:
        dict: A new dictionary containing only entries where the key matches the `filter_key`.
    """
    return {key: dictionary[key] for key in dictionary.keys() if key == filter_key}


def diffs_from_dict_series(
    dict_series: dict[float, dict[str, ProcessResult]],
    dict_thematic: dict[str, BaseGeometry],
    diff_metric: DiffMetric = DiffMetric.CHANGES_AREA,
):
    """
    Calculates a dictionary containing difference metrics for thematic elements based on a distance series.

    This function analyzes the changes in thematic elements (represented by thematic_ids in `dict_thematic`)
    across different distances provided in the `dict_series`. It calculates a difference metric
    for each thematic element at each distance and returns a dictionary summarizing these differences.

    Args:
        dict_series (dict): A dictionary where thematic_ids are distances and values are tuples of two dictionaries.
                             - The first dictionary in the tuple represents thematic element areas for a specific distance.
                             - The second dictionary represents the difference in areas from the original thematic data for a specific distance.
        dict_thematic (dict): A dictionary where thematic_ids are thematic element identifiers and values are GeoJSON geometry objects
                               representing the original thematic data.
       diff_metric (DiffMetric): The metric used to determine the difference between the thematic and reference data.

    Returns:
        dict: A dictionary containing difference metrics for each thematic element (`key`) across different distances.
             The structure is as follows:
             {
                'thematic_key1': {
                    distance1: difference_metric1,
                    distance2: difference_metric2,
                    ...
                },
                'thematic_key2': {
                    distance1: difference_metric1,
                    distance2: difference_metric2,
                    ...
                },
                ...
             }

             - `difference_metric`: This value depends on the chosen calculation for thematic element change.
               The docstring provides examples like area difference, percentage change, and absolute difference.

    Raises:
        KeyError: If a thematic element key is missing from the results in `dict_series`.
    """
    thematic_ids = dict_thematic.keys()

    diffs = {thematic_id: {} for thematic_id in thematic_ids}

    # all the relevant distances used to calculate the series
    for rel_dist, results_dict in dict_series.items():

        for thematic_id in thematic_ids:
            result = results_dict.get(thematic_id, {}).get("result")
            result_diff = results_dict.get(thematic_id, {}).get("result_diff")
            # result_diff_plus = results_dict.get(thematic_id, {}).get("result_diff_plus")
            # result_diff_min = results_dict.get(thematic_id, {}).get("result_diff_min")

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
                diff = (
                    result_diff.area
                )  # equals the symmetrical difference, so equal to result_diff_plus.area + result_diff_min.area
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

    This function retrieves a collection of features from a URL that supports pagination using a `startIndex` parameter.
    It iteratively retrieves features in chunks of the specified `limit` until no more features are available.

    Args:
        ref_url (str): The base URL of the API endpoint.
        limit (int): The maximum number of features to retrieve per request.

    Returns:
        dict: A dictionary representing the complete GeoJSON feature collection. This might be truncated
              if the total number of features exceeds the limitations of the API or server.

    Logs:
        - Debug logs the URL being used for each request during pagination.
    """
    start_index = 0
    collection = {}
    while True:
        url = ref_url + "&startIndex=" + str(start_index)
        logging.debug(url)
        json = requests.get(url).json()
        feature_collection = json
        if (
            "features" not in feature_collection
            or len(feature_collection["features"]) == 0
        ):
            break
        start_index = start_index + limit
        if collection == {}:
            collection = feature_collection
        else:
            collection["features"].extend(feature_collection["features"])
    return collection


def merge_process_results(
    result_dict: dict[str, ProcessResult]
) -> dict[str, ProcessResult]:
    """
    Merges geometries in a dictionary from multiple themes into a single theme.

    Args: result_dict (dict): A dictionary where keys are theme IDs and values are process results

    Returns: dict: A new dictionary with merged geometries, where keys are global
        theme IDs and values are merged geometries.

    """
    grouped_results: dict[str, ProcessResult] = {}

    for id_theme, process_result in result_dict.items():
        id_theme_global = id_theme.split(MULTI_SINGLE_ID_SEPARATOR)[0]
        if id_theme_global not in grouped_results:
            grouped_results[id_theme_global] = process_result
        else:
            for key in process_result:
                geom: BaseGeometry = process_result[key]  # noqa
                if geom.is_empty or geom is None:
                    continue
                existing: BaseGeometry = grouped_results[id_theme_global][key]  # noqa
                grouped_results[id_theme_global][key] = unary_union(  # noqa
                    [existing, geom]
                )

    return grouped_results


def processresult_to_dicts(dict_processresult):
    """
    Transforms a dictionary with all ProcessResults to individual dictionaries of the results
    Args:
        dict_processresult:

    Returns:

    """
    results = {}
    results_diff = {}
    results_diff_plus = {}
    results_diff_min = {}
    results_relevant_intersection = {}
    results_relevant_diff = {}
    for key in dict_processresult:
        processresult = dict_processresult[key]
        results[key] = processresult["result"]
        results_diff[key] = processresult["result_diff"]
        results_diff_plus[key] = processresult["result_diff_plus"]
        results_diff_min[key] = processresult["result_diff_min"]
        results_relevant_intersection[key] = processresult[
            "result_relevant_intersection"
        ]
        results_relevant_diff[key] = processresult["result_relevant_diff"]

    return (
        results,
        results_diff,
        results_diff_plus,
        results_diff_min,
        results_relevant_intersection,
        results_relevant_diff,
    )


def dict_predicted_by_keys(dict_predicted):
    """
    Transforms a dict_predicted into a dictionary with theme_id as keys, and a dictionary with all predicted distances and their resulting geometry as a value.
    Args:
        dict_predicted: a dictionary result of the 'predictor'

    Returns: dictionary with theme_id as keys, and a dictionary with all predicted distances and their resulting geometry as a value.

    """
    dict_predicted_by_keys = {}
    for dist, res in dict_predicted.items():
        for key in dict_predicted[dist]:
            result = {key: res[key]}
            if key not in dict_predicted_by_keys.keys():
                dict_predicted_by_keys[key] = {}
            dict_predicted_by_keys[key][dist] = result
    return dict_predicted_by_keys
