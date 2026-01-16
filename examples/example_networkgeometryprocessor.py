import numpy as np

from brdr.aligner import Aligner
from brdr.be.grb.enums import GRBType
from brdr.be.grb.loader import GRBActualLoader
from brdr.configs import ProcessorConfig
from brdr.enums import AlignerResultType
from brdr.geometry_utils import geom_from_wkt
from brdr.loader import DictLoader
from brdr.processor import NetworkGeometryProcessor
from brdr.viz import plot_difference_by_relevant_distance, show_map

wkt = "MULTIPOLYGON (((161124.34599999999045394 192109.59909999999217689, 161168.89799999998649582 192116.24609999998938292, 161193.2410000000090804 192119.87810000000172295, 161209.3779999999969732 192122.2851000000082422, 161301.39799999998649582 192136.21809999999823049, 161331.3720000000030268 192140.75609999999869615, 161333.79800000000977889 192162.08809999999357387, 161337.36100000000442378 192190.08609999998589046, 161340.79999999998835847 192217.1010999999998603, 161289.46599999998579733 192225.79409999999916181, 161292.97899999999208376 192242.98009999998612329, 161302.7559999999939464 192241.45209999999497086, 161310.06700000001001172 192277.30910000001313165, 161310.77600000001257285 192280.7890999999945052, 161307.41599999999743886 192281.32409999999799766, 161308.2559999999939464 192285.75409999999101274, 161291.07000000000698492 192289.80210000000079162, 161294.72800000000279397 192307.19709999999031425, 161286.91599999999743886 192307.39009999998961575, 161271.16300000000046566 192309.70809999998891726, 161271.12599999998928979 192310.10310000000754371, 161269.97500000000582077 192322.45410000000265427, 161267.09099999998579733 192353.4101000000082422, 161269.3219999999855645 192359.60910000000149012, 161250.23999999999068677 192360.83710000000428408, 161261.48999999999068677 192355.92110000000684522, 161263.70499999998719431 192321.40810000000055879, 161236.57800000000861473 192316.75109999999403954, 161206.75800000000162981 192311.63109999999869615, 161203.84299999999348074 192311.13109999999869615, 161177.42100000000209548 192306.5951000000059139, 161186.53400000001420267 192247.45209999999497086, 161139.78800000000046566 192240.76809999998658895, 161140.12100000001373701 192238.20809999998891726, 161124.34599999999045394 192109.59909999999217689)))"
dict_theme = {"id_wkt": geom_from_wkt(wkt)}
config = ProcessorConfig()
# config.od_strategy = OpenDomainStrategy.EXCLUDE
# processor= DieussaertGeometryProcessor(config=config)
processor = NetworkGeometryProcessor(config=config)
aligner = Aligner(crs="EPSG:31370", processor=processor)
loader = DictLoader(dict_theme)
aligner.load_thematic_data(loader)
loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
aligner.load_reference_data(loader)

series = np.arange(0, 310, 10, dtype=int) / 100
# predict which relevant distances are interesting to propose as resulting geometry
aligner_result = aligner.predict(
    relevant_distances=series,
)
# aligner_result = aligner.process(relevant_distances=[0,1,2,3,4],)
# aligner_result = aligner.process(relevant_distances=[1],)
dict_predictions = aligner_result.get_results(
    aligner=aligner, result_type=AlignerResultType.PREDICTIONS
)

# SHOW results of the predictions
fcs = aligner_result.get_results_as_geojson(add_metadata=False, aligner=aligner)
diffs_dict = aligner.get_difference_metrics_for_thematic_data(
    dict_processresults=aligner_result.results, thematic_data=aligner.thematic_data
)
reference_geometries = {
    key: feat.geometry for key, feat in aligner.reference_data.features.items()
}
if fcs is None or "result" not in fcs:
    print("empty predictions")
else:
    print(fcs["result"])
    for key in dict_predictions:
        plot_difference_by_relevant_distance({key: diffs_dict[key]})
        show_map(
            {key: dict_predictions[key]},
            {key: aligner.thematic_data.features[key].geometry},
            reference_geometries,
        )
