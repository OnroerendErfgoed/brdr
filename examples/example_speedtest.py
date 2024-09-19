import statistics
from datetime import datetime

from brdr.aligner import Aligner
from brdr.loader import GeoJsonFileLoader

# Initiate brdr
aligner = Aligner(relevant_distance=2)
aligner.multi_as_single_modus = True
# Load local thematic data and reference data
# loader = GeoJsonFileLoader(
#     "../tests/testdata/theme.geojson", "theme_identifier"
# )
loader = GeoJsonFileLoader(
    "../tests/testdata/themelayer_not_referenced.geojson", "theme_identifier"
)
aligner.load_thematic_data(loader)
loader = GeoJsonFileLoader("../tests/testdata/reference_leuven.geojson", "capakey")
aligner.load_reference_data(loader)

times = []
for iter in range(1, 3):
    starttime = datetime.now()

    # Example how to use the Aligner
    aligner.predictor()
    fcs = aligner.get_results_as_geojson(formula=True)
    endtime = datetime.now()
    seconds = (endtime - starttime).total_seconds()
    times.append(seconds)
    print(seconds)
print("duration: " + str(times))

print("Min: " + str(min(times)))
print("Max: " + str(max(times)))
print("Mean: " + str(statistics.mean(times)))
print("Median: " + str(statistics.median(times)))
print("Stdv: " + str(statistics.stdev(times)))

# #BEFORE REFACTORING dict_series
# duration: [25.652311, 27.894154, 19.641618, 19.929254, 44.754033, 25.218422, 23.167992, 18.649832, 22.899336, 52.108296]
# Min: 18.649832
# Max: 52.108296
# Mean: 27.9915248
# Median: 24.193207
# Stdv: 11.28891821173264

# #AFTER refactoring
# duration: [21.313991, 16.558168, 16.590126, 18.111118, 16.872433, 17.928071, 18.32295, 17.87116, 19.516652, 16.729241]
# Min: 16.558168
# Max: 21.313991
# Mean: 17.981391
# Median: 17.8996155
# Stdv: 1.504459449440969
