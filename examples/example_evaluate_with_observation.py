import json
from datetime import date

import numpy as np

from brdr.aligner import Aligner
from brdr.be.grb.enums import GRBType
from brdr.be.grb.loader import GRBActualLoader, GRBFiscalParcelLoader
from brdr.be.grb.utils import get_affected_ids_by_grb_change
from brdr.constants import (
    EVALUATION_FIELD_NAME,
    PREDICTION_SCORE,
    RELEVANT_DISTANCE_FIELD_NAME,
)
from brdr.enums import AlignerResultType
from brdr.loader import DictLoader, GeoJsonFileLoader

# Press the green button in the gutter to run the script.
if __name__ == "__main__":
    """
    EXAMPLE of the 'evaluate()-function of 'brdr': This function evaluates thematic objects with a former brdr_observation and compares them with an actual observation; and adds evaluation-properties to the result
    """
    # initiate a base Aligner, to align thematic objects on an older version of the parcels (year 2022)
    base_aligner = Aligner()
    # Load thematic data
    loader = GeoJsonFileLoader("input/themelayer.geojson", "theme_identifier")
    base_aligner.load_thematic_data(loader)
    base_year = "2022"
    name_observation = "brdr_metadata"
    # Load reference data
    base_aligner.load_reference_data(
        GRBFiscalParcelLoader(year=base_year, aligner=base_aligner)
    )
    relevant_distance = 5
    # Align the thematic object on the parcelborders of 2022, to simulate a base-situation
    base_process_result = base_aligner.process(relevant_distances=[relevant_distance])

    # Collect the base-situation (base-geometries and the brdr_observation from that moment
    thematic_dict_observation = {}
    thematic_dict_result = {}
    base_results = base_process_result.results
    for key in base_results:
        thematic_dict_result[key] = base_results[key][relevant_distance]["result"]
        thematic_dict_observation[key] = {
            name_observation: json.dumps(
                base_results[key][relevant_distance]["metadata"]
            )
        }
        print(key + ": " + thematic_dict_result[key].wkt)
        print(key + ": " + str(thematic_dict_observation[key]))

    # (OPTIONAL) Check for changes in the period 2022-now of the reference-parcels (GRB/Flanders-specific function)
    affected = get_affected_ids_by_grb_change(
        thematic_geometries=thematic_dict_result,
        grb_type=GRBType.ADP,
        date_start=date(2022, 1, 1),
        date_end=date.today(),
        one_by_one=False,
        border_distance=relevant_distance,
    )
    if len(affected) == 0:
        print("No affected dicts")
        exit()
    print("Affected_IDs: " + str(affected))

    # Start an aligner to align thematic objects on the actual parcels
    actual_aligner = Aligner()
    # Load the thematic objects (aligned on 2022) and also give the brdr_observation from 2022 as property
    actual_aligner.load_thematic_data(
        DictLoader(
            data_dict=thematic_dict_result,
            data_dict_properties=thematic_dict_observation,
        )
    )
    # Load reference data; the actual parcels
    actual_aligner.load_reference_data(
        GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=actual_aligner)
    )
    # Use the EVALUATE-function
    relevant_distances = np.arange(0, 310, 10, dtype=int) / 100
    aligner_result = actual_aligner.evaluate(
        relevant_distances=relevant_distances,
        thematic_ids=affected,
        metadata_field=name_observation,
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
