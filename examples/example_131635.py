from brdr.aligner import Aligner
from brdr.enums import GRBType
from brdr.grb import GRBActualLoader
from brdr.loader import DictLoader
from brdr.utils import get_oe_dict_by_ids
from examples import print_formula
from examples import show_map

if __name__ == "__main__":
    # EXAMPLE for a thematic Polygon (aanduid_id 131635)

    # Initiate brdr
    aligner = Aligner()
    # Load thematic data & reference data
    dict_theme = get_oe_dict_by_ids([131635])
    loader = DictLoader(dict_theme)
    aligner.load_thematic_data(loader)
    loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
    aligner.load_reference_data(loader)

    # RESULTS
    rel_dist = 2
    dict_results_by_distance = {rel_dist: aligner.process_dict_thematic(rel_dist, 4)}
    # put resulting tuple in a dictionary
    aligner.export_results("output/", formula=True)

    show_map(dict_results_by_distance, aligner.dict_thematic, aligner.dict_reference)
    print_formula(dict_results_by_distance, aligner)
