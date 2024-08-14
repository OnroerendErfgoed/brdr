import os
import unittest

import numpy as np
from shapely import Point
from shapely import from_wkt
from shapely.geometry import Polygon
from shapely.geometry import shape

from brdr.aligner import Aligner
from brdr.enums import OpenbaarDomeinStrategy
from brdr.geometry_utils import buffer_neg_pos
from brdr.geometry_utils import grid_bounds
from brdr.loader import GeoJsonLoader
from brdr.typings import FeatureCollection
from brdr.typings import ProcessResult


class TestAligner(unittest.TestCase):
    def setUp(self):
        # Create a sample geometry for testing
        self.sample_aligner = Aligner()
        self.sample_geom = Polygon([(0, 0), (0, 10), (10, 10), (10, 0)])

    def test_buffer_neg_pos(self):
        # Check if the result of the _buffer_neg_pos gives an equal geometry
        out = buffer_neg_pos(self.sample_geom, 3)
        self.assertEqual(self.sample_geom, out)

    def test_grid_bounds_1(self):
        # Test _grid_bounds function
        delta = 1.0
        grid_partitions = grid_bounds(self.sample_geom, delta)

        # Check if the result is a list of Polygon objects
        self.assertIsInstance(grid_partitions, list)
        for partition in grid_partitions:
            self.assertIsInstance(partition, Polygon)

        # Add more specific tests based on your requirements

    def test_grid_bounds_2(self):
        # Test _grid_bounds function
        delta = 2.0
        grid_partitions = grid_bounds(self.sample_geom, delta)

        # Check if the result is a list of Polygon objects
        self.assertIsInstance(grid_partitions, list)
        for partition in grid_partitions:
            self.assertIsInstance(partition, Polygon)

    def test_export_results(self):
        aligner = Aligner()
        aligner.load_thematic_data_dict(
            {"theme_id_1": from_wkt("POLYGON ((0 0, 0 9, 5 10, 10 0, 0 0))")}
        )
        aligner.load_reference_data_dict(
            {"ref_id_1": from_wkt("POLYGON ((0 1, 0 10,8 10,10 1,0 1))")}
        )
        aligner.process_dict_thematic()
        path = "./tmp/"
        aligner.export_results(path=path)
        filenames = [f"{k}.geojson" for k in ProcessResult.__annotations__]
        for file_name in os.listdir(path):
            os.remove(path + file_name)
            assert file_name in filenames
        os.rmdir(path)

    def test_get_formula_full_intersection(self):
        # Test when intersection equals reference geometry
        key = "a"
        ref_dict = {key: self.sample_geom}
        self.sample_aligner.load_reference_data_dict(ref_dict)
        res = self.sample_aligner.get_formula(self.sample_geom, with_geom=True)
        result = res[key]
        self.assertTrue(result["full"])
        self.assertEqual(result["percentage"], 100)

    def test_get_formula_partial_intersection(self):
        # Test when intersection is partial
        key = "a"
        ref_dict = {key: self.sample_geom.buffer(0.5)}
        self.sample_aligner.load_reference_data_dict(ref_dict)
        res = self.sample_aligner.get_formula(self.sample_geom, with_geom=True)
        result = res[key]
        self.assertFalse(result["full"])
        self.assertGreater(result["percentage"], 0)
        self.assertLess(result["percentage"], 100)

    def test_process_geometry(self):
        # Test if processed geometry is equal to reference geometry
        key_ref = "a"
        ref_dict = {key_ref: self.sample_geom}
        self.sample_aligner.load_reference_data_dict(ref_dict)
        process_result = self.sample_aligner.process_geometry(
            self.sample_geom.buffer(0.5)
        )
        self.assertTrue(
            from_wkt(process_result["result"].wkt).equals(
                from_wkt(self.sample_geom.wkt)
            )
        )
        self.assertFalse(process_result["result_diff"].is_empty)

    def test_predictor(self):
        ##Load thematic data & reference data
        thematic_dict = {"theme_id": from_wkt("POLYGON ((0 0, 0 9, 5 10, 10 0, 0 0))")}
        # ADD A REFERENCE POLYGON TO REFERENCE DICTIONARY
        reference_dict = {"ref_id": from_wkt("POLYGON ((0 1, 0 10,8 10,10 1,0 1))")}
        # LOAD THEMATIC DICTIONARY
        self.sample_aligner.load_thematic_data_dict(thematic_dict)
        # LOAD REFERENCE DICTIONARY
        self.sample_aligner.load_reference_data_dict(reference_dict)
        series = np.arange(0, 300, 10, dtype=int) / 100
        # predict which relevant distances are interesting to propose as resulting geometry

        dict_series, dict_predicted, dict_diffs = self.sample_aligner.predictor(
            relevant_distances=series, od_strategy=4, threshold_overlap_percentage=50
        )
        self.assertEqual(len(dict_predicted[0]), len(thematic_dict))

    def test_load_reference_data_grb_actual_adp(self):
        thematic_dict = {
            "theme_id_1": from_wkt(
                "MultiPolygon (((174184.09476602054201066 171899.68933439542888664, 174400.56834639035514556 171832.959863749332726, 174388.65236948925303295 171770.99678386366576888, 174182.10876987033407204 171836.13745758961886168, 174184.88916448061354458 171873.07698598300339654, 174184.09476602054201066 171899.68933439542888664)))"
            )
        }
        self.sample_aligner.load_thematic_data_dict(thematic_dict)
        # LOAD REFERENCE DICTIONARY
        self.sample_aligner.load_reference_data_grb_actual(
            grb_type="adp", partition=1000
        )
        self.assertGreater(len(self.sample_aligner.dict_reference), 0)

    def test_load_reference_data_grb_actual_gbg(self):
        thematic_dict = {
            "theme_id_1": from_wkt(
                "MultiPolygon (((174184.09476602054201066 171899.68933439542888664, 174400.56834639035514556 171832.959863749332726, 174388.65236948925303295 171770.99678386366576888, 174182.10876987033407204 171836.13745758961886168, 174184.88916448061354458 171873.07698598300339654, 174184.09476602054201066 171899.68933439542888664)))"
            )
        }
        self.sample_aligner.load_thematic_data_dict(thematic_dict)
        # LOAD REFERENCE DICTIONARY
        self.sample_aligner.load_reference_data_grb_actual(
            grb_type="gbg", partition=1000
        )
        self.assertGreater(len(self.sample_aligner.dict_reference), 0)

    def test_load_reference_data_grb_actual_knw(self):
        thematic_dict = {
            "theme_id_1": from_wkt(
                "MultiPolygon (((174184.09476602054201066 171899.68933439542888664, 174400.56834639035514556 171832.959863749332726, 174388.65236948925303295 171770.99678386366576888, 174182.10876987033407204 171836.13745758961886168, 174184.88916448061354458 171873.07698598300339654, 174184.09476602054201066 171899.68933439542888664)))"
            )
        }
        self.sample_aligner.load_thematic_data_dict(thematic_dict)
        # LOAD REFERENCE DICTIONARY
        self.sample_aligner.load_reference_data_grb_actual(grb_type="knw", partition=0)
        self.sample_aligner.process_dict_thematic()
        self.assertGreaterEqual(len(self.sample_aligner.dict_reference), 0)

    def test_all_od_strategies(self):
        ##Load thematic data & reference data
        thematic_dict = {
            "theme_id_1": from_wkt("POLYGON ((0 0, 0 9, 5 10, 10 0, 0 0))")
        }
        # ADD A REFERENCE POLYGON TO REFERENCE DICTIONARY
        reference_dict = {"ref_id_1": from_wkt("POLYGON ((0 1, 0 10,8 10,10 1,0 1))")}
        # LOAD THEMATIC DICTIONARY
        self.sample_aligner.load_thematic_data_dict(thematic_dict)
        # LOAD REFERENCE DICTIONARY
        self.sample_aligner.load_reference_data_dict(reference_dict)
        for od_strategy in OpenbaarDomeinStrategy:
            process_result = self.sample_aligner.process_dict_thematic(
                relevant_distance=1,
                od_strategy=od_strategy,
                threshold_overlap_percentage=50,
            )
            self.assertEqual(len(process_result["theme_id_1"]), 6)

    def test_process_interior_ring(self):
        thematic_dict = {
            "theme_id_1": from_wkt(
                "MultiPolygon (((174370.67910432978533208 171012.2469546866195742, 174461.24052877808571793 171002.71417316573206335, 174429.46459037516615354 170869.25523187351063825, 174373.85669817010057159 170894.6759825958579313, 174369.09030740967136808 170896.6619787460367661, 174370.67910432978533208 171012.2469546866195742),(174400.07184735251939856 170963.78864862219779752, 174401.26344504262669943 170926.4519209987774957, 174419.53460962430108339 170922.47992869839072227, 174429.06739114515949041 170958.62505863170372322, 174400.07184735251939856 170963.78864862219779752)))"
            )
        }
        self.sample_aligner.load_thematic_data_dict(thematic_dict)
        # LOAD REFERENCE DICTIONARY
        self.sample_aligner.load_reference_data_grb_actual(
            grb_type="adp", partition=1000
        )
        result_dict = self.sample_aligner.process_dict_thematic()
        self.assertEqual(len(result_dict), len(thematic_dict))

    def test_process_circle(self):
        geometry = Point(0, 0).buffer(3)
        thematic_dict = {"key": geometry}
        self.sample_aligner.load_thematic_data_dict(thematic_dict)
        # LOAD REFERENCE DICTIONARY
        self.sample_aligner.load_reference_data_grb_actual(
            grb_type="adp", partition=1000
        )
        results_dict = self.sample_aligner.process_dict_thematic()
        self.assertEqual(geometry, results_dict["key"]["result"])

    def test__prepare_thematic_data(self):
        aligner = Aligner()
        geojson: FeatureCollection = {
            "type": "FeatureCollection",
            "name": "theme",
            "crs": {
                "type": "name",
                "properties": {"name": "urn:ogc:def:crs:EPSG::31370"},
            },
            "features": [
                {
                    "type": "Feature",
                    "properties": {"fid": 4, "id": 4, "theme_identifier": "4"},
                    "geometry": {
                        "type": "MultiPolygon",
                        "coordinates": [
                            [
                                [
                                    [174647.924166895216331, 170964.980246312363306],
                                    [174693.204879119351972, 170943.531487890402786],
                                    [174696.382472959638108, 170930.82111252922914],
                                    [174678.905706838035258, 170901.428369506524177],
                                    [174660.634542256389977, 170861.708446502889274],
                                    [174641.56897921464406, 170820.399726579111302],
                                    [174593.905071610264713, 170839.46528962085722],
                                    [174614.559431572153699, 170881.568408004706725],
                                    [174628.064205393398879, 170926.054721768799936],
                                    [174647.924166895216331, 170964.980246312363306],
                                ]
                            ]
                        ],
                    },
                }
            ],
        }
        thematic_loader = GeoJsonLoader(_input=geojson, id_property="theme_identifier")
        aligner.dict_thematic = thematic_loader.load_data()
        assert aligner.dict_thematic == {"4": shape(geojson["features"][0]["geometry"])}
        self.assertGreater(len(aligner.dict_thematic), 0)

    def test_get_predictions_as_dict(self):
        self.sample_aligner.load_thematic_data_dict(
            {"theme_id_1": from_wkt("POLYGON ((0 0, 0 9, 5 10, 10 0, 0 0))")}
        )
        self.sample_aligner.load_reference_data_dict(
            {"ref_id_1": from_wkt("POLYGON ((0 1, 0 10,8 10,10 1,0 1))")}
        )
        self.sample_aligner.process_dict_thematic()
        self.sample_aligner.get_predictions_as_dict()

    def test_get_results_as_dict(self):
        self.sample_aligner.load_thematic_data_dict(
            {
                "theme_id_1": from_wkt("POLYGON ((0 0, 0 9, 5 10, 10 0, 0 0))"),
                "theme_id_2": from_wkt("POLYGON ((0 0, 0 9, 5 10, 10 0, 0 0))"),
            }
        )
        self.sample_aligner.load_reference_data_dict(
            {"ref_id_1": from_wkt("POLYGON ((0 1, 0 10,8 10,10 1,0 1))")}
        )
        self.sample_aligner.process_dict_thematic()
        result = self.sample_aligner.get_results_as_dict()
        assert len(result) == 2

    def test_get_reference_as_geojson(self):
        self.sample_aligner.load_thematic_data_dict(
            {"theme_id_1": from_wkt("POLYGON ((0 0, 0 9, 5 10, 10 0, 0 0))")}
        )
        self.sample_aligner.load_reference_data_dict(
            {"ref_id_1": from_wkt("POLYGON ((0 1, 0 10,8 10,10 1,0 1))")}
        )
        self.sample_aligner.process_dict_thematic()
        self.sample_aligner.get_reference_as_geojson()

    def test_fully_aligned_input(self):
        aligned_shape = from_wkt("POLYGON ((0 0, 0 9, 5 10, 10 0, 0 0))")
        self.sample_aligner.load_thematic_data_dict({"theme_id_1": aligned_shape})
        self.sample_aligner.load_reference_data_dict({"ref_id_1": aligned_shape})
        result = self.sample_aligner.process_dict_thematic()
        assert result["theme_id_1"].get("result") == aligned_shape
        assert result["theme_id_1"].get("result_diff") is None
        assert result["theme_id_1"].get("result_diff_min") is None
        assert result["theme_id_1"].get("result_diff_plus") is None
