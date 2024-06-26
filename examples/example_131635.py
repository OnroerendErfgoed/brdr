from brdr.aligner import Aligner
from brdr.utils import get_oe_dict_by_ids
from examples import show_map

if __name__ == "__main__":
    # EXAMPLE for a thematic Polygon (aanduid_id 131635)

    # Initiate brdr
    aligner = Aligner()
    # Load thematic data & reference data
    dict_theme = get_oe_dict_by_ids([131635])
    aligner.load_thematic_data_dict(dict_theme)
    aligner.load_reference_data_grb_actual(grb_type="adp", partition=1000)

    # RESULTS
    rel_dist = 2
    dict_results_by_distance = {}
    # put resulting tuple in a dictionary
    dict_results_by_distance[rel_dist] = aligner.process_dict_thematic(rel_dist, 4)
    aligner.export_results("output/")
    show_map(dict_results_by_distance, aligner.dict_thematic, aligner.dict_reference)
    for key in dict_results_by_distance[rel_dist][0]:
        print(key)
        print(aligner.get_formula(dict_results_by_distance[rel_dist][0][key]))
