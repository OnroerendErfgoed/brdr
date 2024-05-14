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

def check_business_equality(base_formula, actual_formula):
    """
    function that checks if 2 formulas are equal (determined by business-logic)
    """
    if base_formula.keys() != actual_formula.keys():
        return False
    for key in base_formula.keys():
        if base_formula[key]['full'] != actual_formula[key]['full']:
            return False
        # if abs(base_formula[key]['area'] - actual_formula[key]['area'])>1: #area changed by 1 m²
        #     return False
        # if abs(base_formula[key]['area'] - actual_formula[key]['area'])/base_formula[key]['area'] > 0.01: #area changed by 1%
        #     return False
    return True


def get_collection(ref_url, limit):
    """
    function to get a collection of features from a url
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


#PARAMS
#=========
crs = "EPSG:31370"
limit = 10000
# bbox = "172800,170900,173000,171100"
# bbox = "172000,172000,174000,174000"
bbox = "170000,170000,175000,174900"
# bbox = "100000,195000,105000,195900"
# bbox = "150000,210000,155000,214900"
#bbox = "173500,173500,174000,174000" # example "aanduid_id" = 34195
base_year = "2023"
base_correction = 0 # relevant distance that is used to align the original geometries to the reference-polygons of the base-year
excluded_area = 10000 #geometries bigger than this, will be excluded
series = [0, 0.5, 1, 1.5, 2] #series of relevant distance that is used to check if we can auto-align the geometries to the actual reference-polygons to get an 'equal' formula


#BASE
#=====
# Initiate a AutoReferencer to create a themeset that is base-referenced on a specific base_year
base_auto_referencer = AutoReferencer()
# Load the thematic data to evaluate
# base_auto_referencer.load_thematic_data_file("testdata/theme_changetest.geojson",
# 'theme_identifier') base_auto_referencer.load_thematic_data_file(
# "testdata/theme_leuven.geojson", 'aanduid_id')
base_auto_referencer.load_thematic_data_geojson(
    get_oe_geojson_by_bbox(bbox), "aanduid_id"
)
logging.info(
    "Number of OE-thematic features loaded into base-autoreferencer: "
    + str(len(base_auto_referencer.dict_thematic))
)
# Load the Base reference data
ref_url = (
    "https://geo.api.vlaanderen.be/Adpf/ogc/features/collections/Adpf"
    + base_year
    + "/items?"
    "limit=" + str(limit) + "&crs=" + crs + "&bbox-crs=EPSG:31370&bbox=" + bbox
)
collection = get_collection(ref_url,limit)

base_auto_referencer.load_reference_data_geojson(collection, "CAPAKEY")

# SEARCH FOR CHANGED Parcels in specific timespan
#=================================================
version_date = base_year + "-01-01"  # all changes between this date and NOW will be found
base_url = "https://geo.api.vlaanderen.be/GRB/ogc/features/collections/ADP/items?"
adp_url = (base_url + "datetime="+ version_date + "%2F9999-12-31T00:00:00Z"
           + "& limit="+ str(limit)+ "&crs=" + crs + "&bbox-crs=EPSG:31370&bbox=" + bbox
)

collection = get_collection(adp_url, limit)
array_features = []
for feature in collection['features']:
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

for key in thematic_intersections:
    geometry_base_original = base_auto_referencer.dict_thematic[key]
    print(key)
    print (geometry_base_original)

#ACTUAL
# Initiate a AutoReferencer to reference thematic features to the actual borders
#================================================================================
actual_auto_referencer = AutoReferencer()
actual_url = (
    "https://geo.api.vlaanderen.be/GRB/ogc/features/collections/ADP/items?"
    "limit=" + str(limit) + "&crs=" + crs + "&bbox-crs=EPSG:31370&bbox=" + bbox
)
collection = get_collection(actual_url, limit)
actual_auto_referencer.load_reference_data_geojson(collection, "CAPAKEY")


#LOOP AND PROCESS ALL POSSIBLE AFFECTED FEATURES
#=================================================
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
        geometry_base_original, base_correction
    )
    last_version_date = base_auto_referencer.get_last_version_date(geometry_base_result)
    logging.info('key:' + key + '-->last_version_date: ' + last_version_date)
    logging.info("Original formula: " + key)
    base_formula = base_auto_referencer.get_formula(geometry_base_result)

    for i in series:
        geometry_actual_result, b, c, d, e, f = actual_auto_referencer.process_geometry(
            geometry_base_result, i
        )
        logging.info(
            "New formula: " + key + " with relevant distance(m) : " + str(i)
        )
        actual_formula = actual_auto_referencer.get_formula(geometry_actual_result)
        diff = True
        if check_business_equality(base_formula, actual_formula):  # Logic to be determined by business
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
