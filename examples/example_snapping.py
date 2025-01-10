from brdr.aligner import Aligner
from brdr.enums import GRBType, OpenbaarDomeinStrategy
from brdr.geometry_utils import geom_from_wkt
from brdr.grb import GRBActualLoader
from brdr.loader import DictLoader
from examples import show_map, print_brdr_formula

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
    #themelayer id 650
    wkt = "Polygon((176502.06903269811300561 174525.21342084850766696, 176598.26893461190047674 174575.05957464542007074, 176586.83924329539877363 174386.15217649791156873, 176567.47226634246180765 174389.32709075248567387, 176431.90342767190304585 174383.29475366877159104, 176421.42621063179103658 174522.03850659393356182, 176502.06903269811300561 174525.21342084850766696))"
    loader = DictLoader({"1": geom_from_wkt(wkt)}
    )
    aligner.load_thematic_data(loader)
    # Load reference data
    loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
    aligner.load_reference_data(loader)
    aligner.partial_snapping=False
    dict_results = aligner.process(
        relevant_distances=[10],
        od_strategy=OpenbaarDomeinStrategy.SNAP_ALL_SIDE,
        threshold_overlap_percentage=50,
    )
    show_map(dict_results, aligner.dict_thematic, aligner.dict_reference)
    print_brdr_formula(dict_results, aligner)
