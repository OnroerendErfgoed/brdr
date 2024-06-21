from brdr.aligner import Aligner
from brdr.enums import OpenbaarDomeinStrategy
from brdr.utils import diffs_from_dict_series, geojson_tuple_from_series, write_geojson
from examples import plot_series
from examples import show_map

if __name__ == "__main__":
    # Initiate brdr
    aligner = Aligner()
    # Load thematic data
    aligner.load_thematic_data_file(
        "../tests/testdata/themelayer_referenced.geojson", "id_theme"
    )

    # Use GRB adp-parcels as reference polygons adp= administratieve percelen
    aligner.load_reference_data_grb_actual(grb_type="adp", partition=1000)
    #alternative reference poly
    # # Use GRB-gbg (buildings), gbg= gebouw aan de grond
    # x.load_reference_data_grb_actual('gbg')
    ## Use local data
    # x.load_reference_data_file("../tests/testdata/reference_leuven.geojson", 'capakey')

    # Example how to use the Aligner
    rel_dist = 10
    dict_results_by_distance = {}
    dict_results_by_distance[rel_dist] = aligner.process_dict_thematic(relevant_distance=rel_dist, od_strategy=OpenbaarDomeinStrategy.SNAP_FULL_AREA_SINGLE_SIDE)
    aligner.export_results("output/")
    show_map(dict_results_by_distance, aligner.dict_thematic, aligner.dict_reference)

    rel_dist = 6
    dict_results_by_distance = {}
    dict_results_by_distance[rel_dist] = aligner.process_dict_thematic(relevant_distance=rel_dist, od_strategy=OpenbaarDomeinStrategy.SNAP_ALL_SIDE)
    aligner.export_results("output/")
    show_map(dict_results_by_distance, aligner.dict_thematic, aligner.dict_reference)
    # for key in r:
    #     x.get_formula(r[key])

    # Example how to use a series (for histogram)
    series = [0.1, 0.2, 0.3, 0.4, 0.5, 1, 2]
    dict_series = aligner.process_series(series, 2, 50)
    resulting_areas = diffs_from_dict_series(dict_series, aligner.dict_thematic)
    plot_series(series, resulting_areas)

    # Example how to use the Aligner with treshold_overlap_percentage=-1 (original
    # border will be used for cases where relevant zones cannot be used for determination)
    rel_dist = 6
    dict_results_by_distance = {}
    dict_results_by_distance[rel_dist] = aligner.process_dict_thematic(relevant_distance=rel_dist, od_strategy=OpenbaarDomeinStrategy.SNAP_FULL_AREA_ALL_SIDE,treshold_overlap_percentage=-1)
    aligner.export_results("output/")
    show_map(dict_results_by_distance, aligner.dict_thematic, aligner.dict_reference)