import json
from datetime import date

import numpy as np

from brdr.aligner import Aligner
from brdr.constants import (
    EVALUATION_FIELD_NAME,
    PREDICTION_SCORE,
    RELEVANT_DISTANCE_FIELD_NAME,
)
from brdr.enums import GRBType, AlignerResultType
from brdr.grb import GRBActualLoader
from brdr.grb import GRBFiscalParcelLoader
from brdr.grb import get_affected_by_grb_change
from brdr.loader import DictLoader, GeoJsonFileLoader

# Press the green button in the gutter to run the script.
if __name__ == "__main__":
    """
    EXAMPLE of the 'evaluate()-function of 'brdr': This function evaluates thematic objects with a former brdr_formula and compares them with an actual formula; and adds evaluation-properties to the result
    """
    # initiate a base Aligner, to align thematic objects on an older version of the parcels (year 2022)
    base_aligner = Aligner()
    # Load thematic data
    loader = GeoJsonFileLoader("input/themelayer.geojson", "theme_identifier")
    base_aligner.load_thematic_data(loader)
    base_year = "2022"
    name_formula = "base_formula"
    # Load reference data
    base_aligner.load_reference_data(
        GRBFiscalParcelLoader(year=base_year, aligner=base_aligner)
    )
    relevant_distance = 2
    # Align the thematic object on the parcelborders of 2022, to simulate a base-situation
    base_process_result = base_aligner.process(relevant_distance=relevant_distance)

    # Collect the base-situation (base-geometries and the brdr_formula from that moment
    thematic_dict_formula = {}
    thematic_dict_result = {}
    for key in base_process_result:
        thematic_dict_result[key] = base_process_result[key][relevant_distance][
            "result"
        ]
        thematic_dict_formula[key] = {
            name_formula: json.dumps(
                base_aligner.get_brdr_formula(thematic_dict_result[key])
            )
        }
        print(key + ": " + thematic_dict_result[key].wkt)
        print(key + ": " + str(thematic_dict_formula[key]))

    # (OPTIONAL) Check for changes in the period 2022-now of the reference-parcels (GRB/Flanders-specific function)
    affected, unaffected = get_affected_by_grb_change(
        dict_thematic=thematic_dict_result,
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
    actual_aligner = Aligner(relevant_distances=np.arange(0, 210, 10, dtype=int) / 100)
    # Load the thematic objects (aligned on 2022) and also give the brdr_formula from 2022 as property
    actual_aligner.load_thematic_data(
        DictLoader(
            data_dict=thematic_dict_result, data_dict_properties=thematic_dict_formula
        )
    )
    # Load reference data; the actual parcels
    actual_aligner.load_reference_data(
        GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=actual_aligner)
    )
    # Use the EVALUATE-function
    dict_evaluated, prop_dictionary = actual_aligner.evaluate(
        ids_to_evaluate=affected, base_formula_field=name_formula
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
            + " - score "
            + str(feature["properties"][PREDICTION_SCORE])
            + " - distance "
            + str(feature["properties"][RELEVANT_DISTANCE_FIELD_NAME])
        )
