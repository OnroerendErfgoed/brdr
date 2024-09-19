import logging

from brdr.aligner import Aligner
from brdr.constants import EVALUATION_FIELD_NAME, RELEVANT_DISTANCE_FIELD_NAME
from brdr.grb import GRBFiscalParcelLoader
from brdr.grb import update_to_actual_grb
from brdr.oe import OnroerendErfgoedLoader

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
bbox = [172800, 170900, 173000, 171100]
bbox = [172000, 172000, 174000, 174000]
# bbox = "170000,170000,175000,174900"
# bbox = "100000,195000,105000,195900"
# bbox = "150000,210000,155000,214900"
# bbox = "173500,173500,174000,174000" # example "aanduid_id" = 34195
base_year = "2022"
# relevant distance that is used to align the original geometries to the
# reference-polygons of the base-year
base_correction = 2
# geometries bigger than this, will be excluded
excluded_area = 10000
# max_distance_for_actualisation  of relevant distance that is used to check if we can auto-align the geometries
# to the actual reference-polygons to get an 'equal' formula
max_distance_for_actualisation = 2
# BASE
# =====
# Initiate an Aligner to create a themeset that is base-referenced on a specific
# base_year
base_aligner = Aligner()
# Load the thematic data to evaluate
loader = OnroerendErfgoedLoader(bbox=bbox, partition=0)
base_aligner.load_thematic_data(loader)

logging.info(
    "Number of OE-thematic features loaded into base-aligner: "
    + str(len(base_aligner.dict_thematic))
)
base_aligner.load_reference_data(
    GRBFiscalParcelLoader(year=base_year, aligner=base_aligner)
)

# Exclude objects bigger than specified area
keys_to_exclude = []
nr_features = len(base_aligner.dict_thematic)
for key in base_aligner.dict_thematic:
    if base_aligner.dict_thematic[key].area > excluded_area:
        keys_to_exclude.append(key)
        counter_excluded = counter_excluded + 1
        logging.info(
            "geometrie excluded; bigger than " + str(excluded_area) + ": " + key
        )
for x in keys_to_exclude:
    del base_aligner.dict_thematic[x]

# # Align the features to the base-GRB
base_process_result = base_aligner.process(relevant_distance=base_correction)
# get resulting aligned features on Adpfxxxx, with formula
processresults = base_aligner.get_results_as_geojson(formula=True)
if len(processresults) == 0:
    print("empty processresults")
    exit()
featurecollection_base_result = processresults["result"]

# Update Featurecollection to actual version
fcs = update_to_actual_grb(
    featurecollection_base_result,
    base_aligner.name_thematic_id,
    max_distance_for_actualisation=max_distance_for_actualisation,
)


counter_equality = 0
counter_equality_by_alignment = 0
counter_difference = 0
for feature in fcs["result"]["features"]:
    if EVALUATION_FIELD_NAME in feature["properties"].keys():
        ev = feature["properties"][EVALUATION_FIELD_NAME]
        rd = feature["properties"][RELEVANT_DISTANCE_FIELD_NAME]
        if ev.startswith("equal") and rd == 0:
            counter_equality = counter_equality + 1
        elif ev.startswith("equal") and rd > 0:
            counter_equality_by_alignment = counter_equality_by_alignment + 1
        else:
            counter_difference = counter_difference + 1

print(
    "Features: "
    + str(nr_features)
    + "//Equality: "
    + str(counter_equality)
    + "//Equality by alignment: "
    + str(counter_equality_by_alignment)
    + "//Difference: "
    + str(counter_difference)
    + "//Excluded: "
    + str(counter_excluded)
)
