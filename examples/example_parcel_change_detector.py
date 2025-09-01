import os
from datetime import datetime

from brdr.aligner import Aligner
from brdr.constants import (
    EVALUATION_FIELD_NAME,
    RELEVANT_DISTANCE_FIELD_NAME,
    FORMULA_FIELD_NAME,
)
from brdr.grb import GRBFiscalParcelLoader
from brdr.grb import update_to_actual_grb
from brdr.oe import OnroerendErfgoedLoader
from brdr.utils import write_geojson

if __name__ == "__main__":
    """
    # This code shows an example how the aligner can be used inside a flow of
    # parcel change detection:
    # * it can be used to do a first alignment of the original features
    #       (in this case, based on parcel version adpF2022 (1st January 2022)
    # * it can be used to do a new alignment on the actual version of the parcels adp
    # * it can be used to convert the geometries to a formula, to compare and
    #       evaluate if equality is detected after alignement
    """
counter_excluded = 0
# PARAMS
# =========
crs = "EPSG:31370"
limit = 10
bbox = [172800, 170900, 173000, 171100]
bbox = [172000, 172000, 174000, 174000]
base_year = "2022"
# relevant distance that is used to align the original geometries to the
# reference-polygons of the base-year
base_correction = 2
# geometries bigger than this, will be excluded
excluded_area = 1000
# max_distance_for_actualisation  of relevant distance that is used to check if we can auto-align the geometries
# to the actual reference-polygons to get an 'equal' formula
max_distance_for_actualisation = 2
# BASE
# =====
# Initiate an Aligner to create a themeset that is base-referenced on a specific
# base_year
base_aligner = Aligner()
print("start loading OE-objects")
# Load the thematic data to evaluate
loader = OnroerendErfgoedLoader(bbox=bbox, partition=0)
base_aligner.load_thematic_data(loader)

print(
    "Number of OE-thematic features loaded into base-aligner: "
    + str(len(base_aligner.dict_thematic))
)
base_aligner.load_reference_data(
    GRBFiscalParcelLoader(year=base_year, aligner=base_aligner, partition=1000)
)
print("Reference-data loaded")
# Exclude objects bigger than specified area
keys_to_exclude = []
nr_features = len(base_aligner.dict_thematic)
for key in base_aligner.dict_thematic:
    if base_aligner.dict_thematic[key].area > excluded_area:
        keys_to_exclude.append(key)
        counter_excluded = counter_excluded + 1
        print("geometrie excluded; bigger than " + str(excluded_area) + ": " + str(key))
for x in keys_to_exclude:
    del base_aligner.dict_thematic[x]

# # Align the features to the base-GRB
print("Process base objects")
starttime = datetime.now()
base_process_result = base_aligner.process(relevant_distance=base_correction)
# get resulting aligned features on Adpfxxxx, with formula
processresults = base_aligner.get_results_as_geojson(formula=True)
if len(processresults) == 0:
    print("empty processresults")
    exit()
featurecollection_base_result = processresults["result"]

# Update Featurecollection to actual version
print("Actualise base objects")
fcs = update_to_actual_grb(
    featurecollection_base_result,
    base_aligner.name_thematic_id,
    base_formula_field=FORMULA_FIELD_NAME,
    max_distance_for_actualisation=max_distance_for_actualisation,
)

write_geojson(
    os.path.join("output/", "parcel_change_detector_with.geojson"), fcs["result"]
)

counter_equality = 0
counter_equality_by_alignment = 0
counter_difference = 0
counter_no_change = 0
# TODO:  counter_difference collects al the 'TO_CHECK's' but these are multiple  proposals, so clean up the stats
# TODO: Move this as general output from the updater?
for feature in fcs["result"]["features"]:
    if EVALUATION_FIELD_NAME in feature["properties"].keys():
        ev = feature["properties"][EVALUATION_FIELD_NAME]
        print(ev)
        rd = feature["properties"][RELEVANT_DISTANCE_FIELD_NAME]
        if ev.startswith("equal") and rd == 0:
            counter_equality = counter_equality + 1
        elif ev.startswith("equal") and rd > 0:
            counter_equality_by_alignment = counter_equality_by_alignment + 1
        elif ev.startswith("no_change"):
            counter_no_change = counter_no_change + 1
        else:
            counter_difference = counter_difference + 1

print(
    "Features: "
    + str(nr_features)
    + "//Equality: "
    + str(counter_equality)
    + "//Equality by alignment: "
    + str(counter_equality_by_alignment)
    + "//No change: "
    + str(counter_no_change)
    + "//Difference: "
    + str(counter_difference)
    + "//Excluded: "
    + str(counter_excluded)
)
endtime = datetime.now()
seconds = (endtime - starttime).total_seconds()
print("duration: " + str(seconds))
