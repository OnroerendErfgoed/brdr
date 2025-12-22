from datetime import datetime

import numpy as np

from brdr.aligner import Aligner
from brdr.be.grb.loader import GRBFiscalParcelLoader
from brdr.be.oe.loader import OnroerendErfgoedLoader

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
excluded_area = 500
# max_distance_for_actualisation  of relevant distance that is used to check if we can auto-align the geometries
# to the actual reference-polygons to get an 'equal' observation
max_distance_for_actualisation = 2
# BASE
# =====
# Initiate an Aligner to create a themeset that is base-referenced on a specific
# base_year
base_aligner = Aligner(max_workers=5)
print("start loading OE-objects")
# Load the thematic data to evaluate
loader = OnroerendErfgoedLoader(bbox=bbox, partition=0)
# loader = OnroerendErfgoedLoader(objectids= [276,120380])
base_aligner.load_thematic_data(loader)

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

print(
    "Number of OE-thematic features loaded into base-aligner: "
    + str(len(base_aligner.dict_thematic))
)
base_aligner.load_reference_data(
    GRBFiscalParcelLoader(year=base_year, aligner=base_aligner, partition=1000)
)
print("Reference-data loaded")


# # Align the features to the base-GRB
print("Process base objects")
starttime = datetime.now()
series = np.arange(0, 210, 10, dtype=int) / 100
base_process_result = base_aligner.process(relevant_distances=series)
# get resulting aligned features on Adpfxxxx, with observation
processresults = base_process_result.get_results_as_geojson(
    aligner=base_aligner, add_metadata=True
)
if len(processresults) == 0:
    print("empty processresults")
    exit()
featurecollection_base_result = processresults["result"]


endtime = datetime.now()
seconds = (endtime - starttime).total_seconds()
print("duration of actualisation: " + str(seconds))
