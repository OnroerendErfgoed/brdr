import numpy as np

from brdr.aligner import Aligner
from brdr.enums import GRBType, AlignerResultType
from brdr.grb import GRBActualLoader
from brdr.loader import GeoJsonFileLoader

# Press the green button in the gutter to run the script.
if __name__ == "__main__":
    """
    example to use the predictor-function to automatically predict which resulting
    geometries are interesting to look at (based on detection of breakpoints and
    relevant distances of 'no-change')
    """
    # Initiate an Aligner
    aligner = Aligner()
    # Load thematic data & reference data
    loader = GeoJsonFileLoader(
        "../tests/testdata/test_wanted_changes.geojson", "theme_id"
    )
    aligner.load_thematic_data(loader)
    loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
    aligner.load_reference_data(loader)

    series = np.arange(0, 300, 10, dtype=int) / 100
    # predict which relevant distances are interesting to propose as resulting geometry
    dict_series, dict_predictions, diffs = aligner.predictor(
        relevant_distances=series, od_strategy=4, threshold_overlap_percentage=50
    )
    fcs = aligner.get_results_as_geojson(
        resulttype=AlignerResultType.PREDICTIONS, formula=False
    )
    print(fcs["result"])
    # for key in dict_predictions:
    #     show_map(
    #         {key:dict_predictions[key]},
    #         {key: aligner.dict_thematic[key]},
    #         aligner.dict_reference,
    #     )
