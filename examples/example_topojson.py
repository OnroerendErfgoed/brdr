import numpy as np

from brdr.aligner import Aligner
from brdr.be.grb.enums import GRBType
from brdr.be.grb.loader import GRBActualLoader
from brdr.configs import ProcessorConfig
from brdr.enums import AlignerResultType
from brdr.loader import GeoJsonFileLoader
from brdr.processor import TopologyProcessor

processor = TopologyProcessor(config=ProcessorConfig(), feedback=None)
aligner = Aligner(crs="EPSG:31370", processor=processor)
loader = GeoJsonFileLoader(
    path_to_file="input/topo_parcels.geojson", id_property="CAPAKEY"
)
aligner.load_thematic_data(loader)
aligner.load_reference_data(
    GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
)

# PREDICT the 'stable' relevant distances, for a series of relevant distances
relevant_distances = np.arange(0, 1010, 20, dtype=int) / 100
# predict which relevant distances are interesting to propose as resulting geometry
aligner_result = aligner.predict(
    relevant_distances=relevant_distances,
)
process_results_predictions = aligner_result.get_results(
    aligner=aligner, result_type=AlignerResultType.PREDICTIONS
)
diffs_dict = aligner.get_difference_metrics_for_thematic_data(
    dict_processresults=aligner_result.results, thematic_data=aligner.thematic_data
)
# SHOW results of the predictions
fcs = aligner_result.get_results_as_geojson(aligner=aligner)
if fcs is None or "result" not in fcs:
    print("empty predictions")
else:
    print(fcs["result"])
