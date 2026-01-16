from brdr.aligner import Aligner
from brdr.be.grb.enums import GRBType
from brdr.be.grb.loader import GRBActualLoader
from brdr.constants import (
    EVALUATION_FIELD_NAME,
    PREDICTION_SCORE,
    RELEVANT_DISTANCE_FIELD_NAME,
)
from brdr.enums import AlignerResultType
from brdr.loader import GeoJsonFileLoader

# Press the green button in the gutter to run the script.
if __name__ == "__main__":
    """
    EXAMPLE of the 'evaluate()-function of 'brdr': This function evaluates thematic objects with a former brdr_observation and compares them with an actual observation; and adds evaluation-properties to the result
    """

    # Start an aligner to align thematic objects on the actual parcels
    actual_aligner = Aligner()
    loader = GeoJsonFileLoader("input/themelayer.geojson", "theme_identifier")
    actual_aligner.load_thematic_data(loader)
    # Load reference data; the actual parcels
    actual_aligner.load_reference_data(
        GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=actual_aligner)
    )
    # Use the EVALUATE-function
    aligner_result = actual_aligner.evaluate(
        relevant_distances=None, thematic_ids=None, metadata_field=None
    )
    # SHOW the EVALUATED results
    fc = aligner_result.get_results_as_geojson(
        result_type=AlignerResultType.EVALUATED_PREDICTIONS,
        add_metadata=True,
        add_original_attributes=True,
        aligner=actual_aligner,
    )
    print(fc["result"])

    for feature in fc["result"]["features"]:
        print(
            feature["properties"][actual_aligner.thematic_data.id_fieldname]
            + ": "
            + feature["properties"][EVALUATION_FIELD_NAME]
            + " - score "
            + str(feature["properties"][PREDICTION_SCORE])
            + " - distance "
            + str(feature["properties"][RELEVANT_DISTANCE_FIELD_NAME])
        )
