from brdr.aligner import Aligner
from brdr.enums import GRBType
from brdr.grb import GRBActualLoader
from brdr.oe import OnroerendErfgoedLoader
from examples import print_brdr_formula
from examples import show_map

if __name__ == "__main__":
    # TODO
    # EXAMPLE for a thematic Polygon (aanduid_id 131635)

    # Initiate brdr
    aligner = Aligner()
    # Load thematic data & reference data
    loader = OnroerendErfgoedLoader([131635])
    aligner.load_thematic_data(loader)
    loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
    aligner.load_reference_data(loader)
    ref_geojson = aligner.get_input_as_geojson()

    # RESULTS
    rel_dist = 2
    dict_results = aligner.process(relevant_distance=rel_dist, od_strategy=4)
    fcs = aligner.get_results_as_geojson()
    print(fcs["result"])
    # put resulting tuple in a dictionary
    aligner.save_results("output/", formula=True)

    show_map(dict_results, aligner.dict_thematic, aligner.dict_reference)
    print_brdr_formula(dict_results, aligner)
