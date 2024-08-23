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
    get_collection_grb_fiscal_parcels, evaluate_grb_affected,
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
thematic_dict_formula = {}
i = 0
for key in base_process_result:
    i = i + 1
    thematic_dict[key] = base_process_result[key]["result"]
    thematic_dict_formula[key] = base_aligner.get_formula(thematic_dict[key])
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
counter_equality, counter_equality_by_alignment, counter_difference =evaluate_grb_affected(dict_affected, thematic_dict_formula,series,actual_aligner)

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
