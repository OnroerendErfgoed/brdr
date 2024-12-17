import numpy as np

from brdr.aligner import Aligner
from brdr.enums import GRBType, AlignerResultType, OpenbaarDomeinStrategy
from brdr.grb import GRBActualLoader
from brdr.loader import GeoJsonFileLoader
from examples import show_map, plot_series

# Press the green button in the gutter to run the script.
if __name__ == "__main__":
    """
    EXAMPLE to use the predictor-function to automatically predict which resulting
    geometries are interesting to look at (based on detection of breakpoints and
    relevant distances of 'no-change')
    """
    # Initiate an Aligner
    aligner = Aligner(max_workers=-1)
    # Load thematic data & reference data
    loader = GeoJsonFileLoader(
        "../tests/testdata/test_wanted_changes.geojson", "theme_id"
    )
    aligner.load_thematic_data(loader)
    # Load reference data
    loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
    aligner.load_reference_data(loader)

    # PREDICT the 'stable' relevant distances, for a series of relevant distances
    series = np.arange(0, 310, 10, dtype=int) / 100
    # predict which relevant distances are interesting to propose as resulting geometry
    dict_series, dict_predictions, diffs = aligner.predictor(
        relevant_distances=series,
        od_strategy=OpenbaarDomeinStrategy.SNAP_ALL_SIDE,
        threshold_overlap_percentage=50,
    )

    # SHOW results of the predictions
    fcs = aligner.get_results_as_geojson(
        resulttype=AlignerResultType.PREDICTIONS, formula=False
    )
    print(fcs["result"])
    for key in dict_predictions:
        plot_series(series, {key: diffs[key]})
        show_map(
            {key: dict_predictions[key]},
            {key: aligner.dict_thematic[key]},
            aligner.dict_reference,
        )
