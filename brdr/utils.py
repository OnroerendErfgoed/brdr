import os.path

import numpy as np
import requests
from geojson import Feature
from geojson import FeatureCollection
from geojson import dump
from shapely import GeometryCollection
from shapely import make_valid
from shapely import node
from shapely import polygonize
from shapely.geometry import shape


def export_geojson(path_to_file, dictionary, crs, name_id, multi_to_single=False):
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
    Util function that checks a dictionary of (multi-)polygons and creates a new
        dictionary with only single Polygons.

    The keys of the original dictionary are partially kept:

    * Polygons keep their key
    * MultiPolygons with only one Polygon keep their key
    * MultiPolygons get a suffix-key (added '_index') for every Polygon
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
    Util function that creates a new version of the reference_dictionary based on the
    existing reference-dictionary. This function is useful:

    *   when there are overlapping polygons in the reference-dictionary:
        It creates non-overlapping polygons based on all reference-borders.
    *   when multiple overlapping references will be combined, f.e parcels and buildings
    As this action makes new polygons, the original reference_id is lost and new
    reference_id is added.
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
            print("geom excluded: " + str(g))
    return dict_ref


def get_oe_dict_by_ids(aanduidingsobjecten):
    """
    Get a thematic dictionary from a list of aanduid_id's
    """
    dict_thematic = {}
    headers = {"Accept": "application/json"}
    for a in aanduidingsobjecten:
        url = "https://inventaris.onroerenderfgoed.be/aanduidingsobjecten/" + str(a)
        response = requests.get(url, headers=headers).json()
        key = str(response["id"])
        geom = shape(response["locatie"]["contour"])
        dict_thematic[key] = geom
    return dict_thematic


def get_oe_geojson_by_bbox(bbox, limit=1000):
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

def get_breakpoints(x, y):
  """
  Determine the extremes of a graph.

  Parameters:
    x (numpy.ndarray): The x values of the graph.
    y (numpy.ndarray): The y values of the graph.

  Returns:
    list: A list of the extremes of the graph.
  """

  y = numerical_derivative(x, y)
  extremen = []
  zero_streak = []
  start_streak = None
  streak = 0
  write_zero_streak = False
  last_extreme = 0
  for i in range(1, len(x)-1):
    if round(y[i],2) == 0:
        streak = streak + 1
        if start_streak == None:
            start_streak = x[i]
    elif (streak !=0):
        write_zero_streak = True
    if write_zero_streak or len(x) - 2 == i:
        end_streak = x[i]
        center_streak = start_streak + (end_streak - start_streak)/2
        zero_streak.append((start_streak,end_streak,center_streak, streak,last_extreme))
        streak = 0
        start_streak = None
        write_zero_streak = False
        print('end_streak')
    if round(y[i],2) > 0 and y[i - 1] <= y[i] and y[i]>= y[i + 1]:
        last_extreme = y[i]
        extremen.append((x[i], y[i], "maximum"))

  for extremum in extremen:
      print(f"{extremum[0]:.2f}, {extremum[1]:.2f} ({extremum[2]})")
  for st in zero_streak:
      print(f"{st[0]:.2f} - {st[1]:.2f} -{st[2]:.2f} - {st[3]:.2f} - startextreme {st[4]:.2f} ")
  print('breakpoint' + str(st))
  #plt.plot(series, afgeleide, label='afgeleide-' + str(key))
  return extremen,zero_streak

def numerical_derivative(x, y):
  """
  Calculate the numerical derivative of a graph.

  Parameters:
    x (numpy.ndarray): The x values of the graph.
    y (numpy.ndarray): The y values of the graph.

  Returns:
    numpy.ndarray: The derivative of y with respect to x.
  """

  dx = x[1] - x[0]
  dy = np.diff(y)
  derivative = dy / dx
  derivative = np.insert(derivative, 0, 0)
  #derivative = np.append(derivative, 0)

  return derivative
