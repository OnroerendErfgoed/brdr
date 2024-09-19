from brdr.aligner import Aligner
from brdr.enums import OpenbaarDomeinStrategy
from brdr.geometry_utils import geom_from_wkt
from brdr.loader import DictLoader

# CREATE AN ALIGNER
aligner = Aligner(
    crs="EPSG:31370",
)
# ADD A THEMATIC POLYGON TO THEMATIC DICTIONARY and LOAD into Aligner
thematic_dict = {"theme_id_1": geom_from_wkt("POLYGON ((0 0, 0 9, 5 10, 10 0, 0 0))")}
loader = DictLoader(thematic_dict)
aligner.load_thematic_data(loader)
# ADD A REFERENCE POLYGON TO REFERENCE DICTIONARY and LOAD into Aligner
reference_dict = {"ref_id_1": geom_from_wkt("POLYGON ((0 1, 0 10,8 10,10 1,0 1))")}
loader = DictLoader(reference_dict)
aligner.load_reference_data(loader)
# EXECUTE THE ALIGNMENT
relevant_distance = 1
process_result = aligner.process(
    relevant_distance=relevant_distance,
    od_strategy=OpenbaarDomeinStrategy.SNAP_SINGLE_SIDE,
    threshold_overlap_percentage=50,
)
# PRINT RESULTS IN WKT
print("result: " + process_result["theme_id_1"][relevant_distance]["result"].wkt)
print(
    "added area: "
    + process_result["theme_id_1"][relevant_distance]["result_diff_plus"].wkt
)
print(
    "removed area: "
    + process_result["theme_id_1"][relevant_distance]["result_diff_min"].wkt
)
