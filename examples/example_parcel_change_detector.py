import logging
from datetime import date, timedelta

import numpy as np
from shapely import STRtree
from shapely.geometry import shape

from brdr.aligner import Aligner
from brdr.enums import GRBType
from brdr.grb import (
    get_last_version_date,
    get_geoms_affected_by_grb_change,
    get_collection_grb_fiscal_parcels,
)
from brdr.loader import GeoJsonLoader, GRBActualLoader, DictLoader
from brdr.utils import get_collection
from brdr.utils import get_oe_geojson_by_bbox


# This code shows an example how the aligner can be used inside a flow of
# parcel change detection:
# * it can be used to do a first alignment of the original features
#       (in this case, based on parcel version adpF2023 (1st January 2023)
# * it can be used to do a new alignment on the actual version of the parcels adp
# * it can be used to convert the geometries to a formula, to compare and
#       evaluate if equality is detected after alignement


# TODO: research naar aanduid_id 116448 (equality na 0.5m), 120194 (1m)
def check_business_equality(base_formula, actual_formula):
    """
    function that checks if 2 formulas are equal (determined by business-logic)
    """
    # TODO: research and implementation of following ideas
    # ideas:
    # * If result_diff smaller than 0.x --> automatic update
    # * big polygons: If 'outer ring' has same formula (do net check inner side) --> automatic update
    # ** outer ring can be calculated: 1) nageative buffer 2) original - buffered
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


counter_equality = 0
counter_equality_by_alignment = 0
counter_difference = 0
counter_excluded = 0
# PARAMS
# =========
crs = "EPSG:31370"
limit = 10000
bbox = "172800,170900,173000,171100"
bbox = "172000,172000,174000,174000"
# bbox = "170000,170000,175000,174900"
# bbox = "100000,195000,105000,195900"
# bbox = "150000,210000,155000,214900"
# bbox = "173500,173500,174000,174000" # example "aanduid_id" = 34195
base_year = "2023"
base_correction = 0.2  # relevant distance that is used to align the original geometries to the reference-polygons of the base-year
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
collection_fiscal_parcels = get_collection_grb_fiscal_parcels(base_year, bbox=bbox)

loader = GeoJsonLoader(collection_fiscal_parcels, "CAPAKEY")
base_aligner.load_reference_data(loader)


#Exclude objects bigger than specified area
keys_to_exclude = []
for key in base_aligner.dict_thematic:
    if base_aligner.dict_thematic[key].area > excluded_area:
        keys_to_exclude.append(key)
        counter_excluded = counter_excluded + 1
        logging.info(
            "geometrie excluded; bigger than " + str(excluded_area) + ": " + key
        )
for x in keys_to_exclude:
    del base_aligner.dict_thematic[x]


# Align the features to the base-GRB
base_process_result = base_aligner.process_dict_thematic(
    relevant_distance=base_correction
)

thematic_dict = {}
i = 0
for key in base_process_result:
    i = i + 1
    thematic_dict[key] = base_process_result[key]["result"]
    if i > 500:
        break


dict_affected = get_geoms_affected_by_grb_change(
    thematic_dict,
    grb_type=GRBType.ADP,
    date_start=date.today() - timedelta(days=365),
    date_end=date.today(),
    one_by_one=False,
)

logging.info(
    "Number of possible affected OE-thematic during timespan: "
    + str(len(dict_affected))
)

# ACTUAL
# Initiate a Aligner to reference thematic features to the actual borders
# ================================================================================

# Initiate a Aligner to reference thematic features to the actual borders
actual_aligner = Aligner()
loader = DictLoader(dict_affected)
actual_aligner.load_thematic_data(loader)
loader = GRBActualLoader(grb_type=GRBType.ADP, partition=0, aligner=actual_aligner)
actual_aligner.load_reference_data(loader)


# LOOP AND PROCESS ALL POSSIBLE AFFECTED FEATURES
# =================================================


for key in dict_affected:
    geometry_base_original = dict_affected[key]

    last_version_date = get_last_version_date(geometry_base_original)
    logging.info("key:" + key + "-->last_version_date: " + str(last_version_date))
    logging.info("Original formula: " + key)
    base_formula = base_aligner.get_formula(geometry_base_original)
    logging.info(str(base_formula))

    for i in series:
        actual_process_result = actual_aligner.process_geometry(
            geometry_base_original, i
        )
        logging.info("New formula: " + key + " with relevant distance(m) : " + str(i))
        actual_formula = actual_aligner.get_formula(actual_process_result["result"])
        logging.info(str(actual_formula))
        diff = True
        if check_business_equality(
            base_formula, actual_formula
        ):  # Logic to be determined by business
            if i == 0:
                counter_equality = counter_equality + 1
                logging.info(
                    "equality detected for: " + key + " at distance(m): " + str(i)
                )
            else:
                counter_equality_by_alignment = counter_equality_by_alignment + 1
                logging.info(
                    "equality_by_alignment detected for: "
                    + key
                    + " at distance(m): "
                    + str(i)
                )
            diff = False
            break
    if diff:
        counter_difference = counter_difference + 1
        logging.info("Difference detected for: " + key)
logging.info(
    "Features: "
    + str(len(dict_affected))
    + "//Equality: "
    + str(counter_equality)
    + "//Equality by alignment: "
    + str(counter_equality_by_alignment)
    + "//Difference: "
    + str(counter_difference)
    + "//Excluded: "
    + str(counter_excluded)
)
