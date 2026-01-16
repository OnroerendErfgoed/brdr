import numpy as np

from brdr.aligner import Aligner
from brdr.be.grb.enums import GRBType
from brdr.be.grb.loader import GRBActualLoader
from brdr.configs import ProcessorConfig
from brdr.enums import AlignerResultType
from brdr.geometry_utils import geom_from_wkt
from brdr.loader import DictLoader
from brdr.processor import (
    NetworkGeometryProcessor,
    DieussaertGeometryProcessor,
)
from brdr.viz import plot_difference_by_relevant_distance, show_map

wkt = "MULTIPOLYGON (((172942.43716106086503714 172661.59368339509819634, 172953.8069053897052072 172656.50436926694237627, 172959.11278607649728656 172667.11613064052653499, 172961.27845166294719093 172671.23089525476098061, 172950.99154012731742114 172677.40304217615630478, 172946.22707583714509383 172680.8681071144528687, 172937.67269677069270983 172664.73389849544037133, 172942.43716106086503714 172661.59368339509819634)))"
# TODO check with geometry below (big geometry) if actual grb download is working
# wkt ="POLYGON ((174086.46188367175636813 174996.72320061543723568, 174425.62381920358166099 174440.8744729382742662, 174345.54391775856493041 173828.49875600583618507, 173916.8809159058437217 173319.75585270809824578, 173247.9782097181014251 173197.28070932161062956, 172838.15753761713858694 173437.52041365666082129, 172673.28715228918008506 173640.07545848813606426, 172418.91570064029656351 173828.49875600583618507, 172188.09716118115466088 174148.81836178587400354, 172112.72784217409207486 174525.66495682124514133, 172296.44055725380894728 174968.45970598777057603, 172583.78608596828416921 175354.72746589902089909, 172894.68452687244280241 175467.78144440962933004, 173276.24170434573898092 175533.72959854081273079, 173629.5353871913976036 175401.83329027844592929, 174086.46188367175636813 174996.72320061543723568))"
dict_theme = {"id_wkt": geom_from_wkt(wkt)}
config = ProcessorConfig()
# config.od_strategy = OpenDomainStrategy.EXCLUDE
processor = DieussaertGeometryProcessor(config=config)
processor = NetworkGeometryProcessor(config=config)
# processor = SnapGeometryProcessor(config=config)
aligner = Aligner(crs="EPSG:31370", processor=processor)
loader = DictLoader(dict_theme)
aligner.load_thematic_data(loader)
loader = GRBActualLoader(grb_type=GRBType.ADP, partition=-1, aligner=aligner)
aligner.load_reference_data(loader)

series = np.arange(0, 410, 10, dtype=int) / 100
# predict which relevant distances are interesting to propose as resulting geometry
aligner_result = aligner.predict(
    relevant_distances=series,
)
# aligner_result = aligner.process(relevant_distances=[0,1,2,3,4],)
# aligner_result = aligner.process(relevant_distances=[3],) #weird behaviour at 1.8
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
