import numpy as np

from brdr.aligner import Aligner
from brdr.constants import EVALUATION_FIELD_NAME
from brdr.enums import GRBType, AlignerResultType, FullStrategy
from brdr.geometry_utils import geom_from_wkt
from brdr.grb import GRBActualLoader
from brdr.loader import DictLoader

# Press the green button in the gutter to run the script.
if __name__ == "__main__":
    """
    EXAMPLE of the 'evaluate()-function of 'brdr': This function evaluates thematic objects with a former brdr_formula and compares them with an actual formula; and adds evaluation-properties to the result
    """

    # Start an aligner to align thematic objects on the actual parcels
    actual_aligner = Aligner()
    relevant_distances = np.arange(0, 400, 100, dtype=int) / 100
    wkt = "Polygon ((174906.57643806317355484 179830.59888437716290355, 174719.95761857370962389 179820.51138062097015791, 174538.38255096235661767 179691.89570772959268652, 174442.55126527859829366 179555.71440702106337994, 174364.37311116815544665 179432.14248600779683329, 174576.21069004805758595 179376.66121534878038801, 174783.00451704987790436 179394.31434692209586501, 174906.57643806317355484 179830.59888437716290355))"
    thematic_dict = {"id_1": geom_from_wkt(wkt)}
    loader = DictLoader(data_dict=thematic_dict)
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
        relevant_distances=relevant_distances,
    )

    # SHOW the EVALUATED results
    fc = actual_aligner.get_results_as_geojson(
        resulttype=AlignerResultType.EVALUATED_PREDICTIONS,
        formula=True,
        attributes=True,
    )
    print(str(fc["result"]))

    for feature in fc["result"]["features"]:
        print(
            feature["properties"][actual_aligner.name_thematic_id]
            + ": "
            + feature["properties"][EVALUATION_FIELD_NAME]
        )
