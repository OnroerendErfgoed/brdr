import logging

import numpy as np
from shapely import STRtree
from shapely.geometry import shape

from brdr.aligner import Aligner
from brdr.grb import get_last_version_date
from brdr.loader import GeoJsonLoader
from brdr.utils import get_collection
from brdr.utils import get_oe_geojson_by_bbox


# This code shows an example how the aligner can be used inside a flow of
# parcel change detection:
# * it can be used to do a first alignment of the original features
#       (in this case, based on parcel version adpF2023 (1st January 2023)
# * it can be used to do a new alignment on the actual version of the parcels adp
# * it can be used to convert the geometries to a formula, to compare and
#       evaluate if equality is detected after alignement

# TODO convert this example with GRB-features


def check_business_equality(base_formula, actual_formula):
    """
    function that checks if 2 formulas are equal (determined by business-logic)
    """
    if base_formula.keys() != actual_formula.keys():
        return False
    for key in base_formula.keys():
        if base_formula[key]["full"] != actual_formula[key]["full"]:
            return False
        # if abs(base_formula[key]['area'] - actual_formula[key]['area'])>1: #area changed by 1 m²
        #     return False
        # if abs(base_formula[key]['area'] - actual_formula[key]['area'])/base_formula[key]['area'] > 0.01: #area changed by 1%
        #     return False
    return True


# PARAMS
# =========
crs = "EPSG:31370"
limit = 10000
# bbox = "172800,170900,173000,171100"
# bbox = "172000,172000,174000,174000"
bbox = "170000,170000,175000,174900"
# bbox = "100000,195000,105000,195900"
# bbox = "150000,210000,155000,214900"
# bbox = "173500,173500,174000,174000" # example "aanduid_id" = 34195
base_year = "2023"
base_correction = 0  # relevant distance that is used to align the original geometries to the reference-polygons of the base-year
excluded_area = 10000  # geometries bigger than this, will be excluded
series = [
    0,
    0.5,
    1,
    1.5,
    2,
]  # series of relevant distance that is used to check if we can auto-align the geometries to the actual reference-polygons to get an 'equal' formula

# BASE
# =====
# Initiate a Aligner to create a themeset that is base-referenced on a specific
# base_year
base_aligner = Aligner()
# Load the thematic data to evaluate
loader = GeoJsonLoader(get_oe_geojson_by_bbox(bbox), "aanduid_id")
base_aligner.load_thematic_data(loader)

logging.info(
    "Number of OE-thematic features loaded into base-aligner: "
    + str(len(base_aligner.dict_thematic))
)
# Load the Base reference data
ref_url = (
    "https://geo.api.vlaanderen.be/Adpf/ogc/features/collections/Adpf"
    + base_year
    + "/items?"
    "limit=" + str(limit) + "&crs=" + crs + "&bbox-crs=EPSG:31370&bbox=" + bbox
)
collection = get_collection(ref_url, limit)

# base_aligner.load_reference_data_geojson(collection, "CAPAKEY")
loader = GeoJsonLoader(collection, "CAPAKEY")
base_aligner.load_reference_data(loader)

# SEARCH FOR CHANGED Parcels in specific timespan
# =================================================
version_date = (
    base_year + "-01-01"
)  # all changes between this date and NOW will be found
base_url = "https://geo.api.vlaanderen.be/GRB/ogc/features/collections/ADP/items?"
adp_url = (
    base_url
    + "datetime="
    + version_date
    + "%2F9999-12-31T00:00:00Z"
    + "& limit="
    + str(limit)
    + "&crs="
    + crs
    + "&bbox-crs=EPSG:31370&bbox="
    + bbox
)
print(adp_url)

collection = get_collection(adp_url, limit)
array_features = []
for feature in collection["features"]:
    geom = shape(feature["geometry"]).buffer(0)
    array_features.append(geom)
if len(array_features) == 0:
    logging.info("No detected changes")
    exit()
logging.info("Changed parcels in timespan: " + str(len(array_features)))
thematic_tree = STRtree(list(base_aligner.dict_thematic.values()))
thematic_items = np.array(list(base_aligner.dict_thematic.keys()))
arr_indices = thematic_tree.query(array_features, predicate="intersects")
thematic_intersections = thematic_items.take(arr_indices[1])
thematic_intersections = list(set(thematic_intersections))
logging.info("Number of affected features: " + str(len(thematic_intersections)))

for key in thematic_intersections:
    geometry_base_original = base_aligner.dict_thematic[key]
    print(key)
    print(geometry_base_original)

# ACTUAL
# Initiate a Aligner to reference thematic features to the actual borders
# ================================================================================

# Initiate a Aligner to reference thematic features to the actual borders
actual_aligner = Aligner()
actual_url = (
    "https://geo.api.vlaanderen.be/GRB/ogc/features/collections/ADP/items?"
    "limit=" + str(limit) + "&crs=" + crs + "&bbox-crs=EPSG:31370&bbox=" + bbox
)
collection = get_collection(actual_url, limit)
loader = GeoJsonLoader(collection, "CAPAKEY")
actual_aligner.load_reference_data(loader)


# LOOP AND PROCESS ALL POSSIBLE AFFECTED FEATURES
# =================================================
counter_equality = 0
counter_difference = 0
counter_excluded = 0

for key in thematic_intersections:
    geometry_base_original = base_aligner.dict_thematic[key]
    if geometry_base_original.area > excluded_area:
        counter_excluded = counter_excluded + 1
        logging.info(
            "geometrie excluded; bigger than " + str(excluded_area) + ": " + key
        )
        continue
    base_process_result = base_aligner.process_geometry(
        geometry_base_original, base_correction
    )
    last_version_date = get_last_version_date(base_process_result["result"])
    logging.info("key:" + key + "-->last_version_date: " + str(last_version_date))
    logging.info("Original formula: " + key)
    base_formula = base_aligner.get_formula(base_process_result["result"])

    for i in series:
        actual_process_result = actual_aligner.process_geometry(
            base_process_result["result"], i
        )
        logging.info("New formula: " + key + " with relevant distance(m) : " + str(i))
        actual_formula = actual_aligner.get_formula(actual_process_result["result"])
        diff = True
        if check_business_equality(
            base_formula, actual_formula
        ):  # Logic to be determined by business
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
