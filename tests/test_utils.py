import unittest

import shapely
from shapely import is_empty
from shapely.geometry import Polygon

from brdr.constants import MULTI_SINGLE_ID_SEPARATOR
from brdr.oe import get_oe_dict_by_ids
from brdr.typings import ProcessResult
from brdr.utils import diffs_from_dict_series

# from brdr.utils import filter_dict_by_key
from brdr.utils import get_breakpoints_zerostreak
from brdr.utils import get_collection
from brdr.utils import merge_process_results
from brdr.utils import multipolygons_to_singles
from brdr.utils import polygonize_reference_data


class TestUtils(unittest.TestCase):
    def setUp(self):
        # Create a sample geometry for testing
        self.sample_series = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

    def tearDown(self):
        pass

    def test_get_breakpoints_zerostreak(self):
        sample_diff = [0, 5, 5, 3, 2, 2, 2, 6, 7, 7, 6]
        breakpoints, zerostreaks = get_breakpoints_zerostreak(
            self.sample_series, sample_diff
        )
        assert len(breakpoints) != 0
        assert len(zerostreaks) != 0

    def test_get_breakpoints_zerostreak_no_zerostreaks(self):
        breakpoints, zerostreaks = get_breakpoints_zerostreak(
            self.sample_series, self.sample_series
        )
        assert len(breakpoints) != 0
        assert len(zerostreaks) == 0

    def test_multipolygons_to_singles_empty_dict(self):
        data = {}
        result = multipolygons_to_singles(data)
        self.assertEqual(result, {})

    def test_multipolygons_to_singles_with_point(self):
        geometry1 = shapely.geometry.Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
        geometry2 = shapely.geometry.Point(0, 0)
        data = {"test_id1": geometry1, "test_id2": geometry2}
        result = multipolygons_to_singles(data)
        self.assertEqual(result, {"test_id1": geometry1})

    def test_multipolygons_to_singles_single_polygon(self):
        geometry = shapely.geometry.Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
        data = {"test_id": geometry}
        result = multipolygons_to_singles(data)
        self.assertEqual(result, data)

    def test_multipolygons_to_singles_multipolygon_single_poly(self):
        geometry = shapely.geometry.Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
        data = {"test_id": shapely.geometry.MultiPolygon([geometry])}
        result = multipolygons_to_singles(data)
        self.assertEqual(result, {"test_id": geometry})

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

    def test_get_oe_dict_by_ids(self):
        aanduid_id = 131635
        dict_thematic = get_oe_dict_by_ids([aanduid_id])
        self.assertFalse(is_empty(dict_thematic[str(aanduid_id)]))

    # def test_get_oe_dict_by_ids_erfgoedobject(self):
    #     eo_id = 206363
    #     dict_thematic = get_oe_dict_by_ids([eo_id], oetype=OEType.EO)
    #     self.assertFalse(is_empty(dict_thematic[str(eo_id)]))
    #
    # def test_get_oe_dict_by_ids_empty(self):
    #     dict_thematic = get_oe_dict_by_ids([])
    #     self.assertEqual(dict_thematic, {})
    #
    # def test_get_oe_dict_by_ids_not_existing(self):
    #     aanduid_id = -1
    #     dict_thematic = get_oe_dict_by_ids([aanduid_id])
    #     self.assertEqual(dict_thematic, {})

    # def test_filter_dict_by_key_empty_dict(self):
    #     data = {}
    #     result = filter_dict_by_key(data, "key")
    #     self.assertEqual(result, {})
    #
    # def test_filter_dict_by_key_single_match(self):
    #     data = {"key1": "value1", "key2": "value2"}
    #     result = filter_dict_by_key(data, "key1")
    #     self.assertEqual(result, {"key1": "value1"})
    #
    # def test_filter_dict_by_key_no_match(self):
    #     data = {"key1": "value1", "key2": "value2"}
    #     result = filter_dict_by_key(data, "key3")
    #     self.assertEqual(result, {})

    def test_diffs_from_dict_series_complete(self):
        """Tests diffs_from_dict_series with complete data."""
        # Mock data
        dict_thematic = {
            "theme_id1": Polygon([(0, 0), (10, 0), (10, 10), (0, 10)]),
            "theme_id2": Polygon([(5, 5), (15, 5), (15, 15), (5, 15)]),
        }
        dict_series = {
            "theme_id1": {
                10: {
                    "result": Polygon([(0, 0), (8, 0), (8, 8), (0, 8)]),
                    "result_diff": Polygon([(2, 2), (6, 2), (6, 6), (2, 6)]),
                }
            },
            "theme_id2": {
                10: {
                    "result": Polygon([(7, 7), (13, 7), (13, 13), (7, 13)]),
                    "result_diff": Polygon([(9, 9), (11, 9), (11, 11), (9, 11)]),
                }
            },
        }
        expected_diffs = {"theme_id1": {10: 16.0}, "theme_id2": {10: 4.0}}

        assert expected_diffs == diffs_from_dict_series(
            dict_series.copy(), dict_thematic.copy()
        )

    def test_get_collection(self):
        base_year = 2023
        limit = 100
        crs = "EPSG:31370"
        bbox = "173500,173500,174000,174000"
        ref_url = (
            "https://geo.api.vlaanderen.be/Adpf/ogc/features/collections/Adpf"
            + str(base_year)
            + "/items?"
            "limit=" + str(limit) + "&crs=" + crs + "&bbox-crs=EPSG:31370&bbox=" + bbox
        )
        collection = get_collection(ref_url, limit)
        self.assertTrue("features" in collection.keys())

    def test_merge_process_results(self):
        key_1 = "key" + MULTI_SINGLE_ID_SEPARATOR + "1"
        key_2 = "key" + MULTI_SINGLE_ID_SEPARATOR + "2"
        key_3 = "key_3"
        process_result_1 = ProcessResult()
        process_result_1["result"] = Polygon([(0, 0), (10, 0), (10, 10), (0, 10)])
        process_result_2 = ProcessResult()
        process_result_2["result"] = Polygon([(0, 0), (8, 0), (8, 8), (0, 8)])
        process_result_3 = ProcessResult()
        process_result_3["result"] = Polygon([(0, 0), (8, 0), (8, 8), (0, 8)])
        testdict = {
            key_1: {0: process_result_1},
            key_2: {0: process_result_2},
            key_3: {0: process_result_3},
        }
        merged_testdict = merge_process_results(testdict)
        assert len(merged_testdict.keys()) == 2
