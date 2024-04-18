import os.path

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


def get_oe_geojson_by_ids(aanduidingsobjecten):
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
