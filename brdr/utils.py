import os.path
import numpy as np
import requests
import logging
from geojson import Feature
from geojson import FeatureCollection
from geojson import dump
from shapely import GeometryCollection
from shapely import make_valid
from shapely import node
from shapely import polygonize
from shapely.geometry import shape

def export_geojson(path_to_file, dictionary, crs, name_id, multi_to_single=False):
    """Exports a dictionary of GeoJSON geometries to a GeoJSON file.

    Args:
        path_to_file (str): The path to the output GeoJSON file.
        dictionary (dict): A dictionary where keys are feature identifiers
            and values are GeoJSON geometry objects.
        crs (str): The Coordinate Reference System (CRS) of the geometries.
        name_id (str): The name of the property to use for the feature identifier
            in the output GeoJSON.
        multi_to_single (bool, optional): If True, multipolygon geometries
            will be converted to single polygons before export. Defaults to False.

    Raises:
        OSError: If there is an error creating the directory containing the output file.
    """
    features = []
    if multi_to_single:
        dictionary = multipolygons_to_singles(dictionary)
    for key in dictionary:
        geom = dictionary[key]
        area = geom.area
        perimeter = geom.length
        if area != 0:
            shape_index = perimeter / area
        else:
            shape_index = -1

        features.append(
            Feature(
                geometry=geom,
                properties={
                    name_id: key,
                    "area": area,
                    "perimeter": perimeter,
                    "shape_index": shape_index,
                },
            )
        )
    crs_geojson = {"type": "name", "properties": {"name": crs}}
    feature_collection = FeatureCollection(features, crs=crs_geojson)
    parent = os.path.dirname(path_to_file)
    os.makedirs(parent, exist_ok=True)
    with open(path_to_file, "w") as f:
        dump(feature_collection, f)


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
                new_key = str(key) + "_" + str(i)
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
        dict_ref (dict): A dictionary where keys are identifiers and values are GeoJSON geometries (assumed to be Polygons or MultiPolygons).

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


def get_oe_dict_by_ids(objectids, oetype='aanduidingsobjecten'):
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
        if 'id' in response.keys():
            key = str(response["id"])
            geom = shape(response["locatie"]["contour"])
            dict_thematic[key] = geom
        else:
            logging.warning('object met id ' + str(a) +' werd niet gevonden in ' + oetype)
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
        if round(derivative[i],2) == 0:
            streak = streak + 1
            if start_streak == None:
                start_streak = x[i-1]
        elif (streak !=0):
            write_zero_streak = True
        if start_streak != None and (write_zero_streak or len(x)-1 == i):
            end_streak = x[i-1]
            if derivative[i] == 0:
                end_streak = x[i]
            center_streak = start_streak + (end_streak - start_streak)/2
            zero_streaks.append((start_streak,end_streak,center_streak, streak,last_extreme))
            streak = 0
            start_streak = None
            write_zero_streak = False
            logging.debug('end_streak')
        if i < len(x)-1 and round(derivative[i],2) > 0 and derivative[i - 1] <= derivative[i] and derivative[i]>= derivative[i + 1]:
            last_extreme = derivative[i]
            extremes.append((x[i], derivative[i], "maximum"))
        if i < len(x)-1 and round(derivative[i],2) < 0 and derivative[i - 1] >= derivative[i] and derivative[i]<= derivative[i + 1]:
            last_extreme = derivative[i]
            extremes.append((x[i], derivative[i], "minimum"))
    for extremum in extremes:
        logging.info(f"breakpoints: relevant_distance:{extremum[0]:.2f}, extreme:{extremum[1]:.2f} ({extremum[2]})")
    for st in zero_streaks:
        logging.info(f"zero_streaks: [{st[0]:.2f} - {st[1]:.2f}] - center:{st[2]:.2f} - counter:{st[3]:.2f} - min/max-extreme:{st[4]:.2f} ")
    #plt.plot(series, afgeleide, label='afgeleide-' + str(key))
    return extremes,zero_streaks

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
  #derivative = np.append(derivative, 0)

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
def filter_resulting_series_by_key(resulting_series,filter_key):
    """
    Filters a dictionary of result tuples based on a specific key.

    This function iterates through a dictionary `resulting_series` where values are tuples containing six dictionaries.
    The function creates a new dictionary `filtered_resulting_series` with the same structure. It iterates through each distance key (`dist`) in the original dictionary and performs the following:

    1. **Filters each inner dictionary:** It uses the helper function `_filter_dict_by_key` to filter each of the six dictionaries within the result tuple at the current distance key. The filtering is based on the provided `filter_key`.
    2. **Creates a new filtered tuple:** A new tuple is created with the filtered dictionaries.
    3. **Adds the filtered tuple to the new dictionary:** The new filtered tuple is added to the `filtered_resulting_series` dictionary using the original distance key (`dist`).

    **Important:** This function assumes the structure of the `resulting_series` dictionary and the existence of the helper function `_filter_dict_by_key`.

    Args:
        resulting_series (dict): A dictionary where keys are distances and values are tuples containing six dictionaries.
        filter_key (str): The key value to filter the inner dictionaries by.

    Returns:
        dict: A new dictionary (`filtered_resulting_series`) with the same structure as the original one, but containing filtered inner dictionaries based on the `filter_key`.
    """
    filtered_resulting_series ={}
    for dist in resulting_series:
        result_tuple = resulting_series[dist]
        filtered_tuple = (
            _filter_dict_by_key(result_tuple[0], filter_key),
            _filter_dict_by_key(result_tuple[1], filter_key),
            _filter_dict_by_key(result_tuple[2], filter_key),
            _filter_dict_by_key(result_tuple[3], filter_key),
            _filter_dict_by_key(result_tuple[4], filter_key),
            _filter_dict_by_key(result_tuple[5], filter_key)
        )
        filtered_resulting_series[dist] = filtered_tuple

    return filtered_resulting_series


def diffs_from_dict_series(dict_series, dict_thematic):
    """
    Calculates a dictionary containing difference metrics for thematic elements based on a distance series.

    This function analyzes the changes in thematic elements (represented by keys in `dict_thematic`)
    across different distances provided in the `dict_series`. It calculates a difference metric
    for each thematic element at each distance and returns a dictionary summarizing these differences.

    Args:
        dict_series (dict): A dictionary where keys are distances and values are tuples of two dictionaries.
                             - The first dictionary in the tuple represents thematic element areas for a specific distance.
                             - The second dictionary represents the difference in areas from the original thematic data for a specific distance.
        dict_thematic (dict): A dictionary where keys are thematic element identifiers and values are GeoJSON geometry objects
                               representing the original thematic data.

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
    keys = dict_thematic.keys()
    diffs = {}
    for key in keys:
        diffs[key]={}
    for rel_dist in dict_series: #all the relevant distances used to calculate the series
        results = dict_series[rel_dist][0]
        results_diff =dict_series[rel_dist][1]
        for key in keys:
            if key not in results.keys() or key not in results_diff.keys():
                raise KeyError("No results calculated for theme_id " + str(key))

            #calculate the diffs you want to have
            #diff = results_diff[key].area * 100 / results[key].area #percentage of change
            diff = results[key].area - dict_thematic[key].area #difference (m²) between area of resulting geometry and original geometry
            #diff = abs(results[key].area - dict_thematic[key].area) #absolute difference (m²) between area of resulting geometry and original geometry
            #diff = abs(results[key].area - dict_thematic[key].area)*100/dict_thematic[key].area #absolute difference (%) between area of resulting geometry and original geometry
            #TODO: determine a good diff-value for determination
            diffs[key][rel_dist] = diff
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
        if "features" not in feature_collection or len(feature_collection["features"]) == 0:
            break
        start_index = start_index + limit
        if collection == {}:
            collection = feature_collection
        else:
            collection["features"].extend(feature_collection["features"])
    return collection
