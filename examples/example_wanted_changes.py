import numpy as np

from brdr.aligner import Aligner
from brdr.utils import (
    get_breakpoints_zerostreak,
    diffs_from_dict_series,
    write_geojson,
    get_series_geojson_dict,
)
from examples import plot_series, show_map, print_formula

# Press the green button in the gutter to run the script.
if __name__ == "__main__":
    """
    example to test if we can seperate wanted deviations from unwanted deviations.
    By calculating results with a series of relevant distances we get a graphic.
    There is a difference in the graphic for wanted and unwanted deviations:
    * unwanted deviations: will evoluate more gradually
    * wanted deviations: have a harder breakpoint
    The type of breakpoint, breakpoint-distance and length between breakpoints can gave an indication of wanted and unwanted deviations
    """
    # TODO: future possibilities to detect the 'wanted' breakpoint automatically? (AI-machine learning?)
    ##Initiate a Aligner
    aligner = Aligner()
    ##Load thematic data & reference data
    aligner.load_thematic_data_file(
        "../tests/testdata/test_wanted_changes.geojson", "theme_id"
    )
    aligner.load_reference_data_grb_actual(
        grb_type="adp", partition=1000
    )  # gebruik de actuele adp-percelen adp= administratieve percelen

    # Example how to use the Aligner
    rel_dist = 2
    dict_results_by_distance = {}
    dict_results_by_distance[rel_dist] = aligner.process_dict_thematic(rel_dist, 4)
    aligner.export_results("output/")
    show_map(dict_results_by_distance, aligner.dict_thematic, aligner.dict_reference)
    # Possibility to get the descriptive formula of a thematic feature
    print_formula(dict_results_by_distance, aligner)

    # Example how to use a series (for histogram)
    series = np.arange(0, 500, 10, dtype=int) / 100
    dict_series = aligner.process_series(series, 4, 50)
    resulting_areas = diffs_from_dict_series(dict_series, aligner.dict_thematic)
    # TODO
    # fc = get_series_geojson_dict(dict_series, aligner.CRS, aligner.name_thematic_id)
    # write_geojson("output/series.geojson", fc[0])
    # write_geojson("output/series_diff.geojson", fc[1])
    # write_geojson("output/series_relevant_difference.geojson", fc[5])
    # plot_series(series, resulting_areas)
    # for key in resulting_areas:
    #     if len(resulting_areas[key]) == len(series):
    #         lst_diffs = list(resulting_areas[key].values())
    #         extremes, zero_streak = get_breakpoints_zerostreak(series, lst_diffs)
    #         print(str(key))
    #         for extremum in extremes:
    #             print(f"{extremum[0]:.2f}, {extremum[1]:.2f} ({extremum[2]})")
    #         for st in zero_streak:
    #             print(
    #                 f"{st[0]:.2f} - {st[1]:.2f} -{st[2]:.2f} - {st[3]:.2f} - startextreme {st[4]:.2f} "
    #             )
    #             dict_results_by_distance = {}
    #             dict_results_by_distance[st[0]] = aligner.process_dict_thematic(
    #                 st[0], 4
    #             )
    #
    #             dict_results_by_distance = filter_resulting_series_by_key(
    #                 dict_results_by_distance, key
    #             )
    #             show_map(
    #                 dict_results_by_distance,
    #                 {key: aligner.dict_thematic[key]},
    #                 aligner.dict_reference,
    #             )
