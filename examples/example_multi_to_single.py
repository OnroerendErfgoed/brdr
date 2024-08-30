from brdr.aligner import Aligner
from brdr.utils import get_oe_dict_by_ids
from brdr.utils import multipolygons_to_singles
from examples import show_map, print_formula

# example to Change a dictionary form multipolygon to single before executing the
# aligner. Can be used on the thematic dictionary as the reference dictionary


if __name__ == "__main__":
    # EXAMPLE for a thematic MultiPolygon
    dict_theme = get_oe_dict_by_ids([110082])

    # WITHOUT MULTI_TO_SINGLE
    # Initiate brdr
    aligner = Aligner()
    # Load thematic data & reference data
    # Get a specific feature of OE that exists out of a Multipolygon

    aligner.load_thematic_data_dict(dict_theme)
    aligner.load_reference_data_grb_actual(grb_type=GRBType.GBG, partition=1000)

    rel_dist = 2
    dict_results_by_distance = {}
    dict_results_by_distance[rel_dist] = aligner.process_dict_thematic(rel_dist, 4)
    aligner.export_results("output/")
    show_map(dict_results_by_distance, aligner.dict_thematic, aligner.dict_reference)

    print_formula(dict_results_by_distance, aligner)

    # WITH MULTI_TO_SINGLE
    # Initiate brdr
    aligner = Aligner()
    # Load thematic data & reference data
    # Get a specific feature of OE that exists out of a Multipolygon
    dict_theme = multipolygons_to_singles(dict_theme)
    aligner.load_thematic_data_dict(dict_theme)
    aligner.load_reference_data_grb_actual(grb_type=GRBType.GBG, partition=1000)

    rel_dist = 5
    dict_results_by_distance = {}
    dict_results_by_distance[rel_dist] = aligner.process_dict_thematic(rel_dist, 4)
    aligner.export_results("output/")
    show_map(dict_results_by_distance, aligner.dict_thematic, aligner.dict_reference)

    print_formula(dict_results_by_distance, aligner)
