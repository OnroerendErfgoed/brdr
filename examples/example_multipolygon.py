# Initiate brdr
import tests
from brdr.aligner import Aligner
from brdr.enums import OpenbaarDomeinStrategy
from brdr.utils import multipolygons_to_singles, write_geojson, dict_series_by_keys
from examples import show_map

aligner0 = Aligner()

# Load thematic data

aligner0.load_thematic_data_file(
    "../tests/testdata/multipolygon.geojson", "theme_identifier"
)
aligner0.dict_thematic = multipolygons_to_singles(aligner0.dict_thematic)
aligner0.load_thematic_data_dict(
    aligner0.dict_thematic,
)
# gebruik de actuele adp-percelen adp= administratieve percelen
aligner = Aligner()
aligner.load_thematic_data_dict(
    aligner0.dict_thematic,
)
aligner.load_reference_data_grb_actual(grb_type="adp", partition=1000)

# Example how to use the Aligner
# rel_dist = 2
# dict_results_by_distance = {}
# dict_results_by_distance[aligner.relevant_distance] = aligner.process_dict_thematic(
#     relevant_distance=rel_dist,
#     od_strategy=OpenbaarDomeinStrategy.SNAP_FULL_AREA_ALL_SIDE,
# )
# aligner.export_results("output/")
# show_map(dict_results_by_distance, aligner.dict_thematic, aligner.dict_reference)
dict_series, dict_predicted, diffs = aligner.predictor()
fcs = aligner.get_predictions_as_geojson(series_dict=dict_predicted, formula=True)
aligner.export_results("output/")
write_geojson("output/predicted.geojson", fcs["result"])
write_geojson("output/predicted_diff.geojson", fcs["result_diff"])
