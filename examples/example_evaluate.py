import numpy as np

from brdr.aligner import Aligner
from brdr.be.grb.enums import GRBType
from brdr.be.grb.loader import GRBActualLoader
from brdr.constants import (
    EVALUATION_FIELD_NAME,
    RELEVANT_DISTANCE_FIELD_NAME,
    PREDICTION_SCORE,
)
from brdr.enums import AlignerResultType, FullReferenceStrategy
from brdr.geometry_utils import geom_from_wkt
from brdr.loader import DictLoader

# Press the green button in the gutter to run the script.
if __name__ == "__main__":
    """
    EXAMPLE of the 'evaluate()-function of 'brdr': This function evaluates thematic objects with a former brdr_observation and compares them with an actual observation; and adds evaluation-properties to the result
    """

    # Start an aligner to align thematic objects on the actual parcels
    aligner = Aligner(crs="EPSG:31370")
    # Define the range of relevant distances to search for predictions
    relevant_distances = np.arange(0, 400, 10, dtype=int) / 100
    # Define the input geometry
    input_geometry_wkt = "POLYGON ((174099.53572337731020525 179375.43192003236617893, 174116.23475504826637916 179372.81056040961993858, 174109.5357249012158718 179324.16977629828033969, 174094.68135370552772656 179324.75230065890355036, 174099.53572337731020525 179375.43192003236617893))"
    thematic_dict = {"id_1": geom_from_wkt(input_geometry_wkt)}
    loader = DictLoader(data_dict=thematic_dict)

    # Load the input geometry into the Aligner
    aligner.load_thematic_data(loader)

    # Load reference data (in this case we use a on-the-fly reference loader to load the actual parcels of Flanders (GRB)
    aligner.load_reference_data(
        GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
    )

    # Use the EVALUATE-function
    aligner_result = aligner.evaluate(
        relevant_distances=relevant_distances,
        thematic_ids=None,
        metadata_field=None,
        full_reference_strategy=FullReferenceStrategy.PREFER_FULL_REFERENCE,
    )

    # Get the EVALUATED results as GEOJSON
    fc = aligner_result.get_results_as_geojson(
        result_type=AlignerResultType.EVALUATED_PREDICTIONS,
        add_metadata=True,
        add_original_attributes=True,
        aligner=aligner,
    )
    # Loop these results to print the results
    for feature in fc["result"]["features"]:
        print(
            feature["properties"][aligner.thematic_data.id_fieldname]
            + ": "
            + feature["properties"][EVALUATION_FIELD_NAME]
            + " - score "
            + str(feature["properties"][PREDICTION_SCORE])
            + " - distance "
            + str(feature["properties"][RELEVANT_DISTANCE_FIELD_NAME])
        )
