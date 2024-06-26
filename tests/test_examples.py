import unittest

import numpy as np
from shapely import from_wkt

from brdr.aligner import Aligner
from brdr.enums import OpenbaarDomeinStrategy
from brdr.utils import (
    get_oe_dict_by_ids,
    multipolygons_to_singles,
    diffs_from_dict_series,
    filter_resulting_series_by_key,
    get_breakpoints_zerostreak,
    get_oe_geojson_by_bbox,
)


class TestExamples(unittest.TestCase):

    def test_example_131635(self):
        # EXAMPLE
        aligner = Aligner()
        dict_theme = get_oe_dict_by_ids([131635])
        aligner.load_thematic_data_dict(dict_theme)
        aligner.load_reference_data_grb_actual(grb_type="adp", partition=1000)
        rel_dist = 2
        dict_results_by_distance = {}
        dict_results_by_distance[rel_dist] = aligner.process_dict_thematic(rel_dist, 4)

    def test_example_combined_borders_adp_gbg(self):
        aligner = Aligner()
        dict_theme = get_oe_dict_by_ids([131635])
        aligner.load_thematic_data_dict(dict_theme)
        dict_adp, name_reference_id_adp = aligner.get_reference_data_dict_grb_actual(
            grb_type="adp", partition=1000
        )
        dict_gbg, name_reference_id_gbg = aligner.get_reference_data_dict_grb_actual(
            grb_type="gbg", partition=1000
        )
        dict_adp_gbg = dict_adp
        dict_adp_gbg.update(dict_gbg)  # combine 2 dictionaries
        # make a polygonized version of the reference data with non-overlapping polygons
        dict_ref = dict_adp_gbg
        aligner.load_reference_data_dict(dict_ref)
        rel_dist = 2
        dict_results_by_distance = {}
        dict_results_by_distance[rel_dist] = aligner.process_dict_thematic(rel_dist, 4)
        results = dict_results_by_distance[rel_dist][0]
        for key in results:
            aligner.get_formula(results[key])

    def test_example_multi_to_single(self):
        aligner = Aligner()
        # Load thematic data & reference data
        # Get a specific feature of OE that exists out of a Multipolygon
        dict_theme = get_oe_dict_by_ids([110082])
        dict_theme = multipolygons_to_singles(dict_theme)
        aligner.load_thematic_data_dict(dict_theme)
        aligner.load_reference_data_grb_actual(grb_type="gbg", partition=1000)

        rel_dist = 5
        dict_results_by_distance = {}
        dict_results_by_distance[rel_dist] = aligner.process_dict_thematic(rel_dist, 4)

        results = dict_results_by_distance[rel_dist][0]

        for key in results:
            print(key)
            print(aligner.get_formula(results[key]))

    # def test_example_parcel_change_detector(self):
    #     aligner_x = Aligner()
    #     dict_theme = get_oe_dict_by_ids([131635])
    #     # Load thematic data & reference data (parcels)
    #     aligner_x.load_thematic_data_dict(dict_theme)
    #     aligner_x.load_reference_data_grb_actual(
    #         grb_type="adp", partition=1000
    #     )  # gebruik de actuele adp-percelen adp= administratieve percelen
    #
    #     aligner_y = Aligner()
    #     # Load thematic data & reference data (buildings)
    #     aligner_y.load_thematic_data_dict(dict_theme, 'aanduid_id')
    #     aligner_y.load_reference_data_grb_actual(
    #         grb_type="gbg", partition=1000
    #     )  # gebruik de actuele adp-percelen adp= administratieve percelen
    #
    #     # Example how to use a series (for histogram)
    #     series = np.arange(0, 300, 10, dtype=int)/100
    #     x_dict_series = aligner_x.process_series(series, 4, 50)
    #     x_resulting_areas = diffs_from_dict_series(x_dict_series, aligner_x.dict_thematic)
    #     y_dict_series = aligner_y.process_series(series, 4, 50)
    #     y_resulting_areas = diffs_from_dict_series(y_dict_series, aligner_y.dict_thematic)
    #     # plot_diffs(series,x_resulting_areas)
    #     # plot_diffs(series,y_resulting_areas)
    #     # Make a 1-by-1 comparison for each thematic feature compared to the 2 references (
    #     # parcels and buildings)
    #     for key in x_resulting_areas:
    #         dict_diff = {"x" + str(key): x_resulting_areas[key], "y" + str(key): y_resulting_areas[key]}

    def test_example_wanted_changes(self):
        aligner = Aligner()
        ##Load thematic data & reference data
        dict_theme = get_oe_dict_by_ids([131635])
        aligner.load_thematic_data_dict(dict_theme)
        aligner.load_reference_data_grb_actual(
            grb_type="adp", partition=1000
        )  # gebruik de actuele adp-percelen adp= administratieve percelen
        # Example how to use the Aligner
        rel_dist = 2
        dict_results_by_distance = {}
        dict_results_by_distance[rel_dist] = aligner.process_dict_thematic(rel_dist, 4)

        # Example how to use a series (for histogram)
        series = np.arange(0, 300, 10, dtype=int) / 100
        dict_series = aligner.process_series(series, 4, 50)
        resulting_areas = diffs_from_dict_series(dict_series, aligner.dict_thematic)
        for key in resulting_areas:
            if len(resulting_areas[key]) == len(series):
                lst_diffs = list(resulting_areas[key].values())
                extremes, zero_streak = get_breakpoints_zerostreak(series, lst_diffs)
                print(str(key))
                for extremum in extremes:
                    print(f"{extremum[0]:.2f}, {extremum[1]:.2f} ({extremum[2]})")
                for st in zero_streak:
                    print(
                        f"{st[0]:.2f} - {st[1]:.2f} -{st[2]:.2f} - {st[3]:.2f} - startextreme {st[4]:.2f} "
                    )
                    dict_results_by_distance = {}
                    dict_results_by_distance[st[0]] = aligner.process_dict_thematic(
                        st[0], 4
                    )
                    dict_results_by_distance = filter_resulting_series_by_key(
                        dict_results_by_distance, key
                    )

    def test_example_predictor(self):
        aligner = Aligner()
        ##Load thematic data & reference data
        dict_theme = get_oe_dict_by_ids([131635])
        aligner.load_thematic_data_dict(dict_theme)
        aligner.load_reference_data_grb_actual(
            grb_type="adp", partition=1000
        )  # gebruik de actuele adp-percelen adp= administratieve percelen

        series = np.arange(0, 300, 10, dtype=int) / 100
        # predict which relevant distances are interesting to propose as resulting geometry
        dict_predicted, diffs = aligner.predictor(
            relevant_distances=series, od_strategy=4, treshold_overlap_percentage=50
        )
        for key in dict_predicted.keys():
            continue
