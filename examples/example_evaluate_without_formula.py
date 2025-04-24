import numpy as np

from brdr.aligner import Aligner
from brdr.constants import EVALUATION_FIELD_NAME
from brdr.enums import GRBType, AlignerResultType, FullStrategy
from brdr.grb import GRBActualLoader
from brdr.loader import GeoJsonFileLoader

# Press the green button in the gutter to run the script.
if __name__ == "__main__":
    """
    EXAMPLE of the 'evaluate()-function of 'brdr': This function evaluates thematic objects with a former brdr_formula and compares them with an actual formula; and adds evaluation-properties to the result
    """

    # Start an aligner to align thematic objects on the actual parcels
    actual_aligner = Aligner(relevant_distances=np.arange(0, 310, 10, dtype=int) / 100)
    loader = GeoJsonFileLoader("input/themelayer.geojson", "theme_identifier")
    actual_aligner.load_thematic_data(loader)
    # Load reference data; the actual parcels
    actual_aligner.load_reference_data(
        GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=actual_aligner)
    )
    # Use the EVALUATE-function
    dict_evaluated, prop_dictionary = actual_aligner.evaluate(
        ids_to_evaluate=None,
        base_formula_field=None,
        full_strategy=FullStrategy.PREFER_FULL,
    )

    # SHOW the EVALUATED results
    fc = actual_aligner.get_results_as_geojson(
        resulttype=AlignerResultType.EVALUATED_PREDICTIONS,
        formula=True,
        attributes=True,
    )
    print(fc["result"])

    for feature in fc["result"]["features"]:
        print(
            feature["properties"][actual_aligner.name_thematic_id]
            + ": "
            + feature["properties"][EVALUATION_FIELD_NAME]
        )
