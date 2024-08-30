from brdr.aligner import Aligner
from brdr.enums import OpenbaarDomeinStrategy
from brdr.geometry_utils import geom_from_wkt
from brdr.loader import DictLoader

# CREATE AN ALIGNER
aligner = Aligner(
    relevant_distance=1,
    od_strategy=OpenbaarDomeinStrategy.SNAP_SINGLE_SIDE,
    threshold_overlap_percentage=50,
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
process_result = aligner.process_dict_thematic(relevant_distance=1)
# PRINT RESULTS IN WKT
print("result: " + process_result["theme_id_1"]["result"].wkt)
print("added area: " + process_result["theme_id_1"]["result_diff_plus"].wkt)
print("removed area: " + process_result["theme_id_1"]["result_diff_min"].wkt)
# SHOW RESULTING GEOMETRY AND CHANGES
# from examples import show_map
# show_map(
#     {aligner.relevant_distance:(result, result_diff, result_diff_plus, result_diff_min, relevant_intersection, relevant_diff)},
#     thematic_dict,
#     reference_dict)
