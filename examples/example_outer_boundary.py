from brdr.aligner import Aligner
from brdr.enums import OpenbaarDomeinStrategy, GRBType
from brdr.geometry_utils import geom_from_wkt
from brdr.grb import GRBActualLoader
from brdr.loader import DictLoader
from examples import print_brdr_formula, show_map

aligner = Aligner(max_workers=None)
od_strategy = OpenbaarDomeinStrategy.SNAP_OUTER_SIDE
relevant_distance=10
wkt = "Polygon ((174015.08694170592934825 179025.39916784031083807, 174040.71934808720834553 179031.93985084796440788, 174037.1838437587430235 178986.50862022733781487, 174030.64316075111855753 178982.97311589887249283, 174016.85469387014745735 179002.24161448894301429, 174021.62762471355381422 179005.77711881740833633, 174018.62244603436556645 179008.60552228015149012, 174015.08694170592934825 179025.39916784031083807))"
geom = geom_from_wkt(wkt)
thematic_dict = {"theme_id_1": geom}
loader = DictLoader(thematic_dict)
aligner.load_thematic_data(loader)
#loader = GeoJsonFileLoader("../tests/testdata/reference_leuven.geojson", "capakey")
#aligner.load_reference_data(loader)
loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
aligner.load_reference_data(loader)

dict_results = aligner.process(relevant_distance=relevant_distance, od_strategy=od_strategy)

# show results
aligner.save_results("output/")
show_map(dict_results, aligner.dict_thematic, aligner.dict_reference)
print_brdr_formula(dict_results, aligner)


