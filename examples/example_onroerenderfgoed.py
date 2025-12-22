from brdr.aligner import Aligner
from brdr.be.grb.enums import GRBType
from brdr.be.grb.loader import GRBActualLoader
from brdr.be.oe.loader import OnroerendErfgoedLoader, OEType
from brdr.enums import AlignerResultType
from examples import print_brdr_observation
from examples import show_map

if __name__ == "__main__":
    # EXAMPLE for a thematic Polygon from Onroerend Erfgoed (https://inventaris.onroerenderfgoed.be/aanduidingsobjecten/131635)

    # Initiate brdr
    aligner = Aligner(max_workers=None)
    # Load thematic data from Onroerend Erfgoed
    loader = OnroerendErfgoedLoader(objectids=[5914, 121125], oetype=OEType.AO)
    aligner.load_thematic_data(loader)
    # Load reference data: The actual GRB-parcels
    loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
    aligner.load_reference_data(loader)

    # PROCESS
    dict_results = aligner.process(relevant_distance=2)

    # GET/SHOW results
    aligner.save_results("output/", observation=True)
    show_map(dict_results, aligner.dict_thematic, aligner.dict_reference)
    print_brdr_observation(dict_results, aligner)

    dict_predictions_evaluated = aligner.evaluate()

    # SHOW results of the predictions
    fcs = aligner.get_results_as_geojson(
        resulttype=AlignerResultType.EVALUATED_PREDICTIONS, observation=False
    )
    print(fcs["result"])
    for key in dict_predictions_evaluated:
        show_map(
            {key: dict_predictions_evaluated[key]},
            {key: aligner.dict_thematic[key]},
            aligner.dict_reference,
        )
