from brdr.aligner import Aligner
from brdr.enums import GRBType
from brdr.grb import GRBActualLoader
from brdr.oe import OnroerendErfgoedLoader, OEType
from examples import print_brdr_formula
from examples import show_map

if __name__ == "__main__":
    # EXAMPLE for a thematic Polygon from Onroerend Erfgoed (https://inventaris.onroerenderfgoed.be/aanduidingsobjecten/131635)

    # Initiate brdr
    aligner = Aligner()
    # Load thematic data from Onroerend Erfgoed
    loader = OnroerendErfgoedLoader(objectids=[131635], oetype=OEType.AO)
    aligner.load_thematic_data(loader)
    # Load reference data: The actual GRB-parcels
    loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
    aligner.load_reference_data(loader)

    # PROCESS
    dict_results = aligner.process(relevant_distance=2)

    # GET/SHOW results
    aligner.save_results("output/", formula=True)
    show_map(dict_results, aligner.dict_thematic, aligner.dict_reference)
    print_brdr_formula(dict_results, aligner)
