import logging

import numpy as np
import requests
from shapely import STRtree
from shapely.geometry import shape

from brdr.auto_referencer import AutoReferencer
from brdr.utils import get_oe_geojson_by_bbox

# This code shows an example how the auto_referencer can be used inside a flow of
# parcel change detection:
# * it can be used to do a first alignment of the original features
#       (in this case, based on parcel version adpF2023 (1st January 2023)
# * it can be used to do a new alignment on the actual version of the parcels adp
# * it can be used to convert the geometries to a formula, to compare and
#       evaluate if equality is detected after alignement


crs = "EPSG:31370"
limit = 10000
# bbox = ("174593.90507161026471294,170820.39972657911130227,174696.38247295963810757,"
#         "170964.98024631236330606")
# bbox = ("171593.90507161026471294,170820.39972657911130227,174696.38247295963810757,"
#         "172064.98024631236330606")
# bbox = "172800,170900,173000,171100"
# bbox = "172000,172000,174000,174000"
bbox = "170000,170000,175000,174900"
base_year = "2023"
excluded_area = 20000
# Initiate a AutoReferencer to create a themeset that is base-referenced on a specific
# base_year
base_auto_referencer = AutoReferencer()
# Load the thematic data to evaluate
# base_auto_referencer.load_thematic_data_file("testdata/theme_changetest.geojson",
# 'theme_identifier') base_auto_referencer.load_thematic_data_file(
# "testdata/theme_leuven.geojson", 'aanduid_id')
base_auto_referencer.load_thematic_data_geojson(
    get_oe_geojson_by_bbox(bbox), "aanduid_id"
)
logging.info(
    "Number of OE-thematic features loaded: "
    + str(len(base_auto_referencer.dict_thematic))
)
# Load the Base reference data
ref_url = (
    "https://geo.api.vlaanderen.be/Adpf/ogc/features/collections/Adpf"
    + base_year
    + "/items?"
    "limit=" + str(limit) + "&crs=" + crs + "&bbox-crs=EPSG:31370&bbox=" + bbox
)
start_index = 0
collection = {}
while True:
    url = ref_url + "&startIndex=" + str(start_index)
    print(url)
    json = requests.get(url).json()
    feature_collection = json
    if "features" not in feature_collection or len(feature_collection["features"]) == 0:
        break
    start_index = start_index + limit
    if collection == {}:
        collection = feature_collection
    else:
        collection["features"].extend(feature_collection["features"])

base_auto_referencer.load_reference_data_geojson(collection, "CAPAKEY")

# SEARCH FOR CHANGED Parcels in specific timespan
version_date = "2023-01-01"  # all changes between this date and NOW will be found
base_url = "https://geo.api.vlaanderen.be/GRB/ogc/features/collections/ADP/items?"
adp_url = (
    base_url
    + "datetime="
    + version_date
    + "%2F9999-12-31T00:00:00Z"
    + "&limit="
    + str(limit)
    + "&crs="
    + crs
    + "&bbox-crs=EPSG:31370&bbox="
    + bbox
)
start_index = 0
array_features = []
while True:
    url = adp_url + "&startIndex=" + str(start_index)
    feature_collection = requests.get(url).json()
    if "features" not in feature_collection or len(feature_collection["features"]) == 0:
        break
    features = feature_collection["features"]
    start_index = start_index + limit
    for feature in features:
        geom = shape(feature["geometry"]).buffer(0)
        array_features.append(geom)
    if len(array_features) == 0:
        logging.info("No detected changes")
        exit()
logging.info("Changed parcels in timespan: " + str(len(array_features)))
thematic_tree = STRtree(list(base_auto_referencer.dict_thematic.values()))
thematic_items = np.array(list(base_auto_referencer.dict_thematic.keys()))
arr_indices = thematic_tree.query(array_features, predicate="intersects")
thematic_intersections = thematic_items.take(arr_indices[1])
thematic_intersections = list(set(thematic_intersections))
logging.info("Number of affected features: " + str(len(thematic_intersections)))

# Initiate a AutoReferencer to reference thematic features to the actual borders
actual_auto_referencer = AutoReferencer()
actual_url = (
    "https://geo.api.vlaanderen.be/GRB/ogc/features/collections/ADP/items?"
    "limit=" + str(limit) + "&crs=" + crs + "&bbox-crs=EPSG:31370&bbox=" + bbox
)
start_index = 0
collection = {}
while True:
    url = actual_url + "&startIndex=" + str(start_index)
    print(url)
    feature_collection = requests.get(url).json()
    if "features" not in feature_collection or len(feature_collection["features"]) == 0:
        break
    start_index = start_index + limit
    if collection == {}:
        collection = feature_collection
    else:
        collection["features"].extend(feature_collection["features"])

actual_auto_referencer.load_reference_data_geojson(collection, "CAPAKEY")

counter_equality = 0
counter_difference = 0
counter_excluded = 0
for key in thematic_intersections:
    geometry_base_original = base_auto_referencer.dict_thematic[key]
    if geometry_base_original.area > excluded_area:
        counter_excluded = counter_excluded + 1
        logging.info(
            "geometrie excluded; bigger than " + str(excluded_area) + ": " + key
        )
        continue
    geometry_base_result, *_ = base_auto_referencer.process_geometry(
        geometry_base_original, 0.2
    )
    print(base_auto_referencer.get_last_version_date(geometry_base_result))
    logging.info("Originele formule voor: " + key)
    base_formula = base_auto_referencer.get_formula(geometry_base_result)
    for i in [0, 0.1, 0.2, 0.5, 1, 1.5, 2]:
        geometry_actual_result, b, c, d, e, f = actual_auto_referencer.process_geometry(
            geometry_base_result, i
        )
        logging.info(
            "Nieuwe formule voor: " + key + " met relevante afstand(m) : " + str(i)
        )
        actual_formula = actual_auto_referencer.get_formula(geometry_actual_result)
        diff = True
        if base_formula == actual_formula:  # Logic to be determined by business
            counter_equality = counter_equality + 1
            logging.info("equality detected for: " + key + " at distance(m): " + str(i))
            diff = False
            break
    if diff:
        counter_difference = counter_difference + 1
        logging.info("Difference detected for: " + key)
logging.info(
    "Features: "
    + str(len(thematic_intersections))
    + "//Equality: "
    + str(counter_equality)
    + "//Difference: "
    + str(counter_difference)
    + "//Excluded: "
    + str(counter_excluded)
)
