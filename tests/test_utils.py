import unittest

import pytest
from shapely.geometry import Polygon

from brdr.utils import (
    get_collection,
    get_geometry_difference_metrics_from_processresults,
)
from brdr.utils import polygonize_reference_data


class TestUtils(unittest.TestCase):
    def setUp(self):
        # Create a sample geometry for testing
        self.sample_series = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    def tearDown(self):
        pass

    # def test_get_breakpoints_zerostreak(self):
    #     sample_diff = [0, 5, 5, 3, 2, 2, 2, 6, 7, 7, 6]
    #     breakpoints, zerostreaks = determine_stability(
    #         self.sample_series, sample_diff
    #     )
    #     assert len(breakpoints) != 0
    #     assert len(zerostreaks) != 0
    #
    # def test_get_breakpoints_zerostreak_no_zerostreaks(self):
    #     breakpoints, zerostreaks = determine_stability(
    #         self.sample_series, self.sample_series
    #     )
    #     assert len(breakpoints) != 0
    #     assert len(zerostreaks) == 0

    def test_polygonize_reference_data_no_overlap(self):
        """Tests polygonize_reference_data with non-overlapping polygons."""
        data = {
            "ref1": Polygon([(0, 0), (1, 0), (1, 1), (0, 1)]),
            "ref2": Polygon([(2, 2), (3, 2), (3, 3), (2, 3)]),
        }
        result = polygonize_reference_data(data.copy())
        # Assert expected number of features and new keys
        self.assertEqual(len(result), 2)
        self.assertIn("1", result.keys())
        self.assertIn("2", result.keys())

    def test_polygonize_reference_data_simple_overlap(self):
        """Tests polygonize_reference_data with a simple overlapping case."""
        data = {
            "ref1": Polygon([(0, 0), (2, 0), (2, 1), (1, 1), (0, 0.5)]),
            "ref2": Polygon([(1, 0), (3, 0), (3, 2), (1, 2)]),
        }
        result = polygonize_reference_data(data.copy())
        # Assert expected number of features (might be more than original due to
        # splitting)
        self.assertGreaterEqual(len(result), 2)
        # Assert new keys are used
        for key in result.keys():
            self.assertTrue(key.isdigit())

    def test_diffs_from_process_results_complete(self):
        """Tests diffs_from_process_results with complete data."""
        # Mock data
        geom_thematic = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
        dict_processresult = {
            10: {
                "result": Polygon([(0, 0), (8, 0), (8, 8), (0, 8)]),
                "result_diff": Polygon([(2, 2), (6, 2), (6, 6), (2, 6)]),
            }
        }
        expected_diffs = {10: 16.0}

        assert expected_diffs == get_geometry_difference_metrics_from_processresults(
            dict_processresult, geom_thematic, reference_union=None
        )

    @pytest.mark.usefixtures("callback_grb_response")
    def test_get_collection(self):
        base_year = 2023
        limit = 100
        crs = "EPSG:31370"
        bbox = "173500,173500,174000,174000"
        ref_url = (
            "https://geo.api.vlaanderen.be/Adpf/ogc/features/collections/Adpf"
            + str(base_year)
            + "/items?"
        )
        params = {"limit": limit, "crs": crs, "bbox-crs": crs, "bbox": bbox}
        collection = get_collection(url=ref_url, params=params)
        self.assertTrue("features" in collection.keys())
