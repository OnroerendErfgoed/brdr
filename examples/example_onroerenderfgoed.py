import json

from brdr.aligner import Aligner
from brdr.be.grb.enums import GRBType
from brdr.be.grb.loader import GRBActualLoader
from brdr.be.oe.loader import OnroerendErfgoedLoader, OEType
from brdr.enums import AlignerResultType
from brdr.viz import print_observation_of_aligner_results
from brdr.viz import show_map

if __name__ == "__main__":
    # EXAMPLE for a thematic Polygon from Onroerend Erfgoed (https://inventaris.onroerenderfgoed.be/aanduidingsobjecten/131635)

    # Initiate brdr
    aligner = Aligner()
    aligner.log_metadata = True
    aligner.add_observations = True
    # Load thematic data from Onroerend Erfgoed
    loader = OnroerendErfgoedLoader(
        objectids=[
            #'https://id.erfgoed.net/aanduidingsobjecten/5914',
            "https://id.erfgoed.net/aanduidingsobjecten/163287"
        ],
        oetype=OEType.AO,
    )
    aligner.load_thematic_data(loader)
    # Load reference data: The actual GRB-parcels
    loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
    aligner.load_reference_data(loader)

    # PROCESS
    aligner_result = aligner.process(relevant_distances=[2])

    # GET/SHOW results
    aligner_result.save_results(
        aligner=aligner, path="output/", add_original_attributes=True, add_metadata=True
    )
    print_observation_of_aligner_results(
        aligner_result.get_results(aligner=aligner), aligner
    )

    thematic_geometries = {
        key: feat.geometry for key, feat in aligner.thematic_data.features.items()
    }
    reference_geometries = {
        key: feat.geometry for key, feat in aligner.reference_data.features.items()
    }
    show_map(aligner_result.results, thematic_geometries, reference_geometries)

    aligner_result = aligner.evaluate()

    # SHOW results of the predictions
    fcs = aligner_result.get_results_as_geojson(
        result_type=AlignerResultType.EVALUATED_PREDICTIONS,
        aligner=aligner,
        add_metadata=True,
    )
    print(
        json.dumps(
            json.loads(fcs["result"]["features"][0]["properties"]["brdr_metadata"]),
            indent=2,
        )
    )
