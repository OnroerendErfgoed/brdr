from brdr.aligner import Aligner
from shapely.geometry import Polygon
from examples import show_results

aligner = Aligner()
thematic_dict= {"theme_id_1": Polygon([(0, 0), (0, 9), (5, 10), (10, 0)])}
reference_dict = {"ref_id_1": Polygon([(0, 1), (0, 10), (8, 10), (10, 1)])}
aligner.load_thematic_data_dict(thematic_dict)
aligner.load_reference_data_dict(reference_dict)
r, rd, rd_plus, rd_min, sd, si = aligner.process_dict_thematic(relevant_distance=1)
show_results(r, rd)

