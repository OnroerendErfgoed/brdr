from brdr.aligner import Aligner
from shapely import from_wkt
from brdr.enums import OpenbaarDomeinStrategy
from examples import show_results

#CREATE AN ALIGNER
aligner = Aligner(relevant_distance=1, od_strategy=OpenbaarDomeinStrategy.SNAP_SINGLE_SIDE,
                  threshold_overlap_percentage=50, crs='EPSG:31370')
#ADD A THEMATIC POLYGON TO THEMATIC DICTIONARY
thematic_dict= {"theme_id_1": from_wkt('POLYGON ((0 0, 0 9, 5 10, 10 0, 0 0))')}
#ADD A REFERENCE POLYGON TO REFERENCE DICTIONARY
reference_dict = {"ref_id_1": from_wkt('POLYGON ((0 1, 0 10,8 10,10 1,0 1))')}
#LOAD THEMATIC DICTIONARY
aligner.load_thematic_data_dict(thematic_dict)
#LOAD REFERENCE DICTIONARY
aligner.load_reference_data_dict(reference_dict)
#EXECUTE THE ALIGNMENT
r, rd, rd_plus, rd_min, sd, si = aligner.process_dict_thematic(relevant_distance=1)
#SHOW RESULTING GEOMETRY AND CHANGES
show_results(r, rd_plus,rd_min,thematic_dict,reference_dict)

