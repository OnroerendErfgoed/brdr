from brdr.aligner import Aligner
from shapely import from_wkt
from brdr.enums import OpenbaarDomeinStrategy

# CREATE AN ALIGNER
aligner = Aligner(
    relevant_distance=1,
    od_strategy=OpenbaarDomeinStrategy.SNAP_SINGLE_SIDE,
    threshold_overlap_percentage=50,
    crs="EPSG:31370",
)
# ADD A THEMATIC POLYGON TO THEMATIC DICTIONARY and LOAD into Aligner
thematic_dict = {"theme_id_1": from_wkt("POLYGON ((0 0, 0 9, 5 10, 10 0, 0 0))")}
aligner.load_thematic_data_dict(thematic_dict)
# ADD A REFERENCE POLYGON TO REFERENCE DICTIONARY and LOAD into Aligner
reference_dict = {"ref_id_1": from_wkt("POLYGON ((0 1, 0 10,8 10,10 1,0 1))")}
aligner.load_reference_data_dict(reference_dict)
# EXECUTE THE ALIGNMENT
(
    result,
    result_diff,
    result_diff_plus,
    result_diff_min,
    relevant_intersection,
    relevant_diff,
) = aligner.process_dict_thematic(relevant_distance=1)
# PRINT RESULTS IN WKT
print("result: " + result["theme_id_1"].wkt)
print("added area: " + result_diff_plus["theme_id_1"].wkt)
print("removed area: " + result_diff_min["theme_id_1"].wkt)
# SHOW RESULTING GEOMETRY AND CHANGES
# from examples import show_map
# show_map(
#     {aligner.relevant_distance:(result, result_diff, result_diff_plus, result_diff_min, relevant_intersection, relevant_diff)},
#     thematic_dict,
#     reference_dict)
