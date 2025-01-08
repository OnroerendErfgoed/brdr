import numpy as np

from brdr.aligner import Aligner
from brdr.enums import GRBType, AlignerResultType, OpenbaarDomeinStrategy
from brdr.geometry_utils import geom_from_wkt
from brdr.grb import GRBActualLoader
from brdr.loader import DictLoader
from examples import show_map, plot_series, print_brdr_formula

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
    wkt = "Polygon ((173935.19719148115837015 179289.70768110611243173, 173933.29809016973013058 179278.39939602438244037, 173938.6501029564824421 179277.10455422112136148, 173936.57835607128799893 179268.99021225408068858, 173940.72184984170598909 179268.12698438524967059, 173949.35412853004527278 179266.14156028692377731, 173953.49762230046326295 179264.93304127056035213, 173964.11532508712843992 179283.92405438493005931, 173935.19719148115837015 179289.70768110611243173))"

    loader = DictLoader({"1": geom_from_wkt(wkt)}
    )
    aligner.load_thematic_data(loader)
    # Load reference data
    loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
    aligner.load_reference_data(loader)

    # PREDICT the 'stable' relevant distances, for a series of relevant distances
    series = np.arange(0, 1010, 10, dtype=int) / 100
    # predict which relevant distances are interesting to propose as resulting geometry
    # dict_series, dict_predictions, diffs = aligner.predictor(
    #     relevant_distances=series,
    #     od_strategy=OpenbaarDomeinStrategy.SNAP_ALL_SIDE,
    #     threshold_overlap_percentage=50,
    # )
    #
    # # SHOW results of the predictions
    # fcs = aligner.get_results_as_geojson(
    #     resulttype=AlignerResultType.PREDICTIONS, formula=False
    # )
    # print(fcs["result"])
    # for key in dict_predictions:
    #     plot_series(series, {key: diffs[key]})
    #     show_map(
    #         {key: dict_predictions[key]},
    #         {key: aligner.dict_thematic[key]},
    #         aligner.dict_reference,
    #     )

    dict_results = aligner.process(
        relevant_distances=[9.5],
        od_strategy=OpenbaarDomeinStrategy.SNAP_ALL_SIDE,
        threshold_overlap_percentage=50,
    )
    show_map(dict_results, aligner.dict_thematic, aligner.dict_reference)
    print_brdr_formula(dict_results, aligner)
