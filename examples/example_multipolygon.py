# Initiate brdr
from brdr.aligner import Aligner
from brdr.enums import GRBType
from brdr.grb import GRBActualLoader
from brdr.loader import DictLoader, GeoJsonFileLoader
from brdr.utils import multipolygons_to_singles
from brdr.utils import write_geojson

aligner0 = Aligner()

# Load thematic data

aligner0.load_thematic_data(
    GeoJsonFileLoader("../tests/testdata/multipolygon.geojson", "theme_identifier")
)
aligner0.dict_thematic = multipolygons_to_singles(aligner0.dict_thematic)
aligner0.load_thematic_data(
    DictLoader(
        aligner0.dict_thematic,
    )
)
# gebruik de actuele adp-percelen adp= administratieve percelen
aligner = Aligner()
aligner.load_thematic_data(DictLoader(aligner0.dict_thematic))

aligner.load_reference_data(
    GRBActualLoader(aligner=aligner, grb_type=GRBType.ADP, partition=1000)
)

dict_series, dict_predicted, diffs = aligner.predictor()
fcs = aligner.get_series_as_geojson(series_dict=dict_predicted, formula=True)
aligner.export_results("output/")
write_geojson("output/predicted.geojson", fcs["result"])
write_geojson("output/predicted_diff.geojson", fcs["result_diff"])
