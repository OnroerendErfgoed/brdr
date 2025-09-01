import numpy as np

from brdr.aligner import Aligner
from brdr.enums import GRBType, OpenDomainStrategy, AlignerResultType
from brdr.grb import GRBActualLoader
from brdr.loader import GeoJsonFileLoader
from examples import plot_series, show_map

aligner = Aligner(crs="EPSG:31370", preserve_topology=True)
loader = GeoJsonFileLoader(
    path_to_file="input/topo_parcels.geojson", id_property="CAPAKEY"
)
aligner.load_thematic_data(loader)
aligner.load_reference_data(
    GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
)
# relevant_distance = 2
# process_result = aligner.process(
#     relevant_distance=relevant_distance,
# )


# PREDICTIONS
# PREDICT the 'stable' relevant distances, for a series of relevant distances
series = np.arange(0, 210, 20, dtype=int) / 100
# predict which relevant distances are interesting to propose as resulting geometry
dict_series, dict_predictions, diffs = aligner.predictor(
    relevant_distances=series,
    od_strategy=OpenDomainStrategy.SNAP_ALL_SIDE,
    threshold_overlap_percentage=50,
)
# SHOW results of the predictions
fcs = aligner.get_results_as_geojson(
    resulttype=AlignerResultType.PREDICTIONS, formula=False
)
if fcs is None or "result" not in fcs:
    print("empty predictions")
else:
    print(fcs["result"])
    for key in dict_predictions:
        plot_series(series, {key: diffs[key]})
        show_map(
            {key: dict_predictions[key]},
            {key: aligner.dict_thematic[key]},
            aligner.dict_reference,
        )

# aligner = Aligner(crs="EPSG:31370", preserve_topology=True)
# loader = GeoJsonFileLoader(
#     path_to_file="input/one_simple.geojson", id_property="CAPAKEY"
# )
# aligner.load_thematic_data(loader)
#
# aligner.load_reference_data(
#     GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
# )
# relevant_distance = 5
# process_result = aligner.process(
#     relevant_distance=relevant_distance,
# )
# print(process_result)
#
# # PREDICT the 'stable' relevant distances, for a series of relevant distances
# series = np.arange(0, 210, 20, dtype=int) / 100
# # predict which relevant distances are interesting to propose as resulting geometry
# dict_series, dict_predictions, diffs = aligner.predictor(
#     relevant_distances=series,
# )
#
#
# # SHOW results of the predictions
# fcs = aligner.get_results_as_geojson(
#     resulttype=AlignerResultType.PREDICTIONS, formula=True
# )
# if fcs is None or "result" not in fcs:
#     print("empty predictions")
# else:
#     print(fcs["result"])
#     for key in dict_predictions:
#         plot_series(series, {key: diffs[key]})
#         show_map(
#             {key: dict_predictions[key]},
#             {key: aligner.dict_thematic[key]},
#             aligner.dict_reference,
#         )
