from brdr.aligner import Aligner
from brdr.enums import OpenbaarDomeinStrategy, GRBType
from brdr.grb import GRBActualLoader
from brdr.loader import GeoJsonFileLoader
from brdr.utils import diffs_from_dict_series
from examples import plot_series
from examples import show_map

if __name__ == "__main__":
    # Initiate brdr
    aligner = Aligner()
    # Load thematic data
    aligner.load_thematic_data(
        GeoJsonFileLoader("../tests/testdata/themelayer_referenced.geojson", "id_theme")
    )

    # Use GRB adp-parcels as reference polygons adp= administratieve percelen
    aligner.load_reference_data(
        GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
    )

    # Example how to use the Aligner
    rel_dist = 10
    dict_results = aligner.process(
        relevant_distance=rel_dist,
        od_strategy=OpenbaarDomeinStrategy.SNAP_FULL_AREA_SINGLE_SIDE,
    )
    aligner.save_results("output/")
    show_map(dict_results, aligner.dict_thematic, aligner.dict_reference)

    rel_dist = 6
    dict_results = aligner.process(
        relevant_distance=rel_dist, od_strategy=OpenbaarDomeinStrategy.SNAP_ALL_SIDE
    )

    aligner.save_results("output/")
    show_map(dict_results, aligner.dict_thematic, aligner.dict_reference)
    # for key in r:
    #     x.get_formula(r[key])

    # Example how to use a series (for histogram)
    series = [0.1, 0.2, 0.3, 0.4, 0.5, 1, 2]
    dict_series = aligner.process(series, 2, 50)
    resulting_areas = diffs_from_dict_series(dict_series, aligner.dict_thematic)
    plot_series(series, resulting_areas)

    # Example how to use the Aligner with threshold_overlap_percentage=-1 (original
    # border will be used for cases where relevant zones cannot be used for
    # determination)
    rel_dist = 6
    dict_results = aligner.process(
        relevant_distance=rel_dist,
        od_strategy=OpenbaarDomeinStrategy.SNAP_FULL_AREA_ALL_SIDE,
        threshold_overlap_percentage=-1,
    )

    aligner.save_results("output/")
    show_map(dict_results, aligner.dict_thematic, aligner.dict_reference)
