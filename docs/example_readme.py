from shapely import from_wkt
from brdr.aligner import Aligner
from shapely.geometry import Polygon
from examples import show_results

#CREATE AN ALIGNER
aligner = Aligner()
#CREATE AN SAMPLE THEMATIC POLYGON
thematic_dict= {"theme_id_1": from_wkt('POLYGON ((0 0, 0 9, 5 10, 10 0, 0 0))')}
#CREATE AN SAMPLE REFERENCE POLYGON
reference_dict = {"ref_id_1": from_wkt('POLYGON ((0 1, 0 10,8 10,10 1,0 1))')}
#LOAD THEMATIC DATA
aligner.load_thematic_data_dict(thematic_dict)
#LOAD REFERENCE DATA
aligner.load_reference_data_dict(reference_dict)
#ALIGN THEMATIC DATA TO REFERENCE DATA
r, rd, rd_plus, rd_min, sd, si = aligner.process_dict_thematic(relevant_distance=1)
#SHOW RESULTING GEOMETRY (BLUE) AND DIFFERENCE (BLACK)
show_results(r, rd_plus,rd_min)

