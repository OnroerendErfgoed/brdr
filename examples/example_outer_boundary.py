from brdr.aligner import Aligner
from brdr.be.grb.enums import GRBType
from brdr.be.grb.loader import GRBActualLoader
from brdr.configs import ProcessorConfig
from brdr.enums import OpenDomainStrategy
from brdr.geometry_utils import geom_from_wkt
from brdr.loader import DictLoader, GeoJsonFileLoader
from brdr.processor import DieussaertGeometryProcessor
from brdr.viz import print_observation_of_aligner_results, show_map

processor_config = ProcessorConfig()
processor_config.od_strategy = OpenDomainStrategy.SNAP_ALL_SIDE
processor = DieussaertGeometryProcessor(config=processor_config)
aligner = Aligner(processor=processor)

relevant_distance = 2
grb_loader = True
wkt = "Polygon ((174015.08694170592934825 179025.39916784031083807, 174040.71934808720834553 179031.93985084796440788, 174037.1838437587430235 178986.50862022733781487, 174030.64316075111855753 178982.97311589887249283, 174016.85469387014745735 179002.24161448894301429, 174021.62762471355381422 179005.77711881740833633, 174018.62244603436556645 179008.60552228015149012, 174015.08694170592934825 179025.39916784031083807))"
# wkt = "Polygon ((173842.69639981063664891 179402.24934623859007843, 173845.30499411537311971 179403.73997155559482053, 173841.76475898752687499 179409.14348832963150926, 173838.59718018895364366 179417.90091206692159176, 173854.9940586757438723 179424.60872599334106781, 173862.07452893140725791 179426.65833580418257043, 173867.47804570547305048 179414.91966143294121139, 173870.08664001018041745 179409.32981649425346404, 173862.81984158989507705 179409.88880098814843222, 173842.69639981063664891 179402.24934623859007843))"
# GROOT
# wkt = "Polygon ((172594.16427666315576062 177732.38409853109624237, 176756.43182838105713017 177257.73955315977218561, 176464.34287738331477158 172310.48294563542003743, 171316.27511604799656197 172675.59413438258343376, 171316.27511604799656197 172675.59413438258343376, 171772.66410198199446313 174446.38339980645105243, 172594.16427666315576062 177732.38409853109624237))"
# grb_loader =False
# wkt = "Polygon ((171681.61391718150116503 174093.33527241024421528, 171699.81010965688619763 174061.65692696082987823, 171718.85482924254029058 174031.39279336185427383, 171679.72830138093559071 174030.5442662516143173, 171669.92309921802370809 174047.60908924668910913, 171681.61391718150116503 174093.33527241024421528))"
geom = geom_from_wkt(wkt)
thematic_dict = {"theme_id_1": geom}
loader = DictLoader(thematic_dict)
aligner.load_thematic_data(loader)
if not grb_loader:
    loader = GeoJsonFileLoader("../tests/testdata/reference_leuven.geojson", "capakey")
else:
    loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
aligner.load_reference_data(loader)

aligner_result = aligner.process(relevant_distances=[relevant_distance])

# show results
aligner_result.save_results(
    aligner=aligner, path="output/", add_original_attributes=True, add_metadata=True
)
print_observation_of_aligner_results(
    aligner_result.get_results(aligner=aligner), aligner
)

thematic_geometries = {
    key: feat.geometry for key, feat in aligner.thematic_data.features.items()
}
reference_geometries = {
    key: feat.geometry for key, feat in aligner.reference_data.features.items()
}
show_map(aligner_result.results, thematic_geometries, reference_geometries)
