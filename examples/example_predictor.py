import numpy as np

from brdr.aligner import Aligner
from brdr.enums import GRBType, AlignerResultType, OpenDomainStrategy
from brdr.geometry_utils import geom_from_wkt
from brdr.grb import GRBActualLoader
from brdr.loader import GeoJsonFileLoader, DictLoader
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

    wkt = "Polygon ((174906.57643806317355484 179830.59888437716290355, 174719.95761857370962389 179820.51138062097015791, 174538.38255096235661767 179691.89570772959268652, 174442.55126527859829366 179555.71440702106337994, 174364.37311116815544665 179432.14248600779683329, 174576.21069004805758595 179376.66121534878038801, 174783.00451704987790436 179394.31434692209586501, 174906.57643806317355484 179830.59888437716290355))"
    thematic_dict = {"id_1": geom_from_wkt(wkt)}
    loader = DictLoader(data_dict=thematic_dict)

    aligner.load_thematic_data(loader)
    # Load reference data
    loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
    aligner.load_reference_data(loader)

    # PREDICT the 'stable' relevant distances, for a series of relevant distances
    series = np.arange(0, 310, 10, dtype=int) / 100
    # predict which relevant distances are interesting to propose as resulting geometry
    dict_series, dict_predictions, diffs = aligner.predictor(
        relevant_distances=series,
        od_strategy=OpenDomainStrategy.SNAP_ALL_SIDE,
        threshold_overlap_percentage=50,
    )

    # SHOW results of the predictions
    fcs = aligner.get_results_as_geojson(
        resulttype=AlignerResultType.PREDICTIONS, formula=False
    )
    if fcs is None or "result" not in fcs:
        print("empty predictions")
    else:
        print(fcs["result"])
        for key in dict_predictions:
            plot_series(series, {key: diffs[key]})
            show_map(
                {key: dict_predictions[key]},
                {key: aligner.dict_thematic[key]},
                aligner.dict_reference,
            )
