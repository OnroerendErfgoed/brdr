import os
from datetime import datetime

from brdr.aligner import Aligner
from brdr.be.grb.grb import update_featurecollection_to_actual_grb
from brdr.be.grb.loader import GRBFiscalParcelLoader
from brdr.be.oe.loader import OnroerendErfgoedLoader
from brdr.constants import (
    EVALUATION_FIELD_NAME,
    RELEVANT_DISTANCE_FIELD_NAME,
    OBSERVATION_FIELD_NAME,
)
from brdr.utils import write_geojson

if __name__ == "__main__":
    """
    # This code shows an example how the aligner can be used inside a flow of
    # parcel change detection:
    # * it can be used to do a first alignment of the original features
    #       (in this case, based on parcel version adpF2022 (1st January 2022)
    # * it can be used to do a new alignment on the actual version of the parcels adp
    # * it can be used to convert the geometries to a observation, to compare and
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
excluded_area = 10000
# max_distance_for_actualisation  of relevant distance that is used to check if we can auto-align the geometries
# to the actual reference-polygons to get an 'equal' observation
max_distance_for_actualisation = 2
# BASE
# =====
# Initiate an Aligner to create a themeset that is base-referenced on a specific
# base_year
base_aligner = Aligner()
print("start loading OE-objects")
# Load the thematic data to evaluate
loader = OnroerendErfgoedLoader(bbox=bbox, partition=1000)
# loader = OnroerendErfgoedLoader(objectids= [276,120380])
base_aligner.load_thematic_data(loader)

# Exclude objects bigger than specified area
keys_to_exclude = []
nr_features = len(base_aligner.thematic_data.features)
# for key in base_aligner.base_aligner.thematic_data.features.keys():
#     if base_aligner.base_aligner.base_aligner.thematic_data.features.get(key).geometry.area > excluded_area:
#         keys_to_exclude.append(key)
#         counter_excluded = counter_excluded + 1
#         print("geometrie excluded; bigger than " + str(excluded_area) + ": " + str(key))
# for x in keys_to_exclude:
#     del base_aligner.dict_thematic[x]

print(
    "Number of OE-thematic features loaded into base-aligner: "
    + str(len(base_aligner.thematic_data.features))
)
base_aligner.load_reference_data(
    GRBFiscalParcelLoader(year=base_year, aligner=base_aligner, partition=1000)
)
print("Reference-data loaded")


# # Align the features to the base-GRB
print("Process base objects")
starttime = datetime.now()
base_process_result = base_aligner.process(relevant_distances=[base_correction])
# get resulting aligned features on Adpfxxxx, with observation
processresults = base_process_result.get_results_as_geojson(
    aligner=base_aligner, add_metadata=True
)
if len(processresults) == 0:
    print("empty processresults")
    exit()
featurecollection_base_result = processresults["result"]

# Update Featurecollection to actual version
starttime = datetime.now()
print(starttime)
print("Actualise base objects")
fcs = update_featurecollection_to_actual_grb(
    featurecollection_base_result,
    # base_aligner.name_thematic_id,
    base_metadata_field=OBSERVATION_FIELD_NAME,
    max_distance_for_actualisation=max_distance_for_actualisation,
)

write_geojson(
    os.path.join("output/", "parcel_change_detector_with.geojson"), fcs["result"]
)

counter_equality = 0
counter_equality_by_alignment = 0
counter_prediction_unique = 0
counter_difference = 0
counter_no_change = 0

# #270: Move this as general output from the updater?
for feature in fcs["result"]["features"]:
    if EVALUATION_FIELD_NAME in feature["properties"].keys():
        ev = feature["properties"][EVALUATION_FIELD_NAME]
        print(ev)
        rd = feature["properties"][RELEVANT_DISTANCE_FIELD_NAME]
        if ev.startswith("equal") and rd == 0:
            counter_equality = counter_equality + 1
        elif ev.startswith("equal") and rd > 0:
            counter_equality_by_alignment = counter_equality_by_alignment + 1
        elif ev.startswith("prediction_unique"):
            counter_prediction_unique = counter_prediction_unique + 1
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
    + "//Prediction (unique): "
    + str(counter_prediction_unique)
    + "//No change: "
    + str(counter_no_change)
    + "//Difference: "
    + str(counter_difference)
    + "//Excluded: "
    + str(counter_excluded)
)
endtime = datetime.now()
seconds = (endtime - starttime).total_seconds()
print(endtime)
print("duration of actualisation: " + str(seconds))
