import numpy as np
from shapely import from_wkt

from brdr.aligner import Aligner
from brdr.enums import GRBType
from brdr.grb import GRBActualLoader
from brdr.loader import DictLoader

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
    loader = DictLoader(
        {
            "id1": from_wkt(
                "MultiPolygon Z (((138430.4033999964594841 194082.86080000177025795 0, 138422.19659999758005142 194080.36510000005364418 0, 138419.01550000160932541 194079.34930000081658363 0, 138412.59849999845027924 194077.14139999821782112 0, 138403.65579999983310699 194074.06430000066757202 0, 138402.19910000264644623 194077.67480000108480453 0, 138401.83420000225305557 194078.57939999923110008 0, 138400.89329999685287476 194080.91140000149607658 0, 138400.31650000065565109 194080.67880000174045563 0, 138399.27300000190734863 194083.37680000066757202 0, 138405.93310000002384186 194085.95410000160336494 0, 138413.51049999892711639 194088.80620000138878822 0, 138427.25680000334978104 194094.29969999939203262 0, 138430.4033999964594841 194082.86080000177025795 0)))"
            )
        }
    )
    aligner.load_thematic_data(loader)
    loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
    aligner.load_reference_data(loader)

    series = np.arange(0, 800, 10, dtype=int) / 100
    # predict which relevant distances are interesting to propose as resulting geometry
    dict_series, dict_predictions, diffs = aligner.predictor(
        relevant_distances=series, od_strategy=4, threshold_overlap_percentage=50
    )
    fcs = aligner.get_results_as_geojson(formula=False)
    print(fcs["result"])
    # for key in dict_predictions:
    #     show_map(
    #         {key:dict_predictions[key]},
    #         {key: aligner.dict_thematic[key]},
    #         aligner.dict_reference,
    #     )
