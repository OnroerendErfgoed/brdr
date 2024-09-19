import os
import unittest
from datetime import date

import numpy as np
from shapely import Point
from shapely import from_wkt
from shapely.geometry import Polygon
from shapely.geometry import shape

from brdr.aligner import Aligner
from brdr.constants import FORMULA_FIELD_NAME
from brdr.enums import GRBType, AlignerResultType
from brdr.enums import OpenbaarDomeinStrategy
from brdr.geometry_utils import _grid_bounds
from brdr.geometry_utils import buffer_neg_pos
from brdr.grb import (
    GRBActualLoader,
    GRBFiscalParcelLoader,
    get_geoms_affected_by_grb_change,
)
from brdr.loader import GeoJsonLoader, DictLoader
from brdr.typings import FeatureCollection, ProcessResult
from brdr.utils import get_series_geojson_dict


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
        grid_partitions = _grid_bounds(self.sample_geom, delta)

        # Check if the result is a list of Polygon objects
        self.assertIsInstance(grid_partitions, list)
        for partition in grid_partitions:
            self.assertIsInstance(partition, Polygon)

        # Add more specific tests based on your requirements

    def test_grid_bounds_2(self):
        # Test _grid_bounds function
        delta = 2.0
        grid_partitions = _grid_bounds(self.sample_geom, delta)

        # Check if the result is a list of Polygon objects
        self.assertIsInstance(grid_partitions, list)
        for partition in grid_partitions:
            self.assertIsInstance(partition, Polygon)

    def test_export_results(self):
        aligner = Aligner()
        aligner.load_thematic_data(
            DictLoader(
                {"theme_id_1": from_wkt("POLYGON ((0 0, 0 9, 5 10, 10 0, 0 0))")}
            )
        )
        aligner.load_reference_data(
            DictLoader({"ref_id_1": from_wkt("POLYGON ((0 1, 0 10,8 10,10 1,0 1))")})
        )
        aligner.process()
        path = "./tmp/"
        resulttype = AlignerResultType.PROCESSRESULTS
        aligner.save_results(path=path, resulttype=resulttype)
        filenames = [
            resulttype.value + f"_{k}.geojson" for k in ProcessResult.__annotations__
        ]
        for file_name in os.listdir(path):
            os.remove(path + file_name)
            assert file_name in filenames
        os.rmdir(path)

    def test_get_formula_full_intersection(self):
        # Test when intersection equals reference geometry
        key = "a"
        ref_dict = {key: self.sample_geom}
        self.sample_aligner.load_reference_data(DictLoader(ref_dict))
        res = self.sample_aligner.get_brdr_formula(self.sample_geom, with_geom=True)
        self.assertTrue(res["full"])
        result = res["reference_features"][key]
        self.assertTrue(result["full"])
        self.assertEqual(result["percentage"], 100)

    def test_get_formula_partial_intersection(self):
        # Test when intersection is partial
        key = "a"
        ref_dict = {key: self.sample_geom.buffer(0.5)}
        self.sample_aligner.load_reference_data(DictLoader(ref_dict))
        res = self.sample_aligner.get_brdr_formula(self.sample_geom, with_geom=True)
        self.assertFalse(res["full"])
        result = res["reference_features"][key]
        self.assertFalse(result["full"])
        self.assertGreater(result["percentage"], 0)
        self.assertLess(result["percentage"], 100)

    def test_process_geometry(self):
        # Test if processed geometry is equal to reference geometry
        key_ref = "a"
        ref_dict = {key_ref: self.sample_geom}
        self.sample_aligner.load_reference_data(DictLoader(ref_dict))
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
        # Load thematic data & reference data
        thematic_dict = {"theme_id": from_wkt("POLYGON ((0 0, 0 9, 5 10, 10 0, 0 0))")}
        # ADD A REFERENCE POLYGON TO REFERENCE DICTIONARY
        reference_dict = {"ref_id": from_wkt("POLYGON ((0 1, 0 10,8 10,10 1,0 1))")}
        # LOAD THEMATIC DICTIONARY
        self.sample_aligner.load_thematic_data(DictLoader(thematic_dict))
        # LOAD REFERENCE DICTIONARY
        self.sample_aligner.load_reference_data(DictLoader(reference_dict))
        series = np.arange(0, 300, 10, dtype=int) / 100
        # predict which relevant distances are interesting to propose as resulting
        # geometry

        dict_series, dict_predictions, dict_diffs = self.sample_aligner.predictor(
            relevant_distances=series, od_strategy=4, threshold_overlap_percentage=50
        )
        self.assertEqual(len(dict_predictions), len(thematic_dict))

    def test_predictor_double_prediction(self):
        """
        Test if a double prediction is filtered out of the prediction results.
        This testdata has 2 resulting predictions that are the same (at 0.0 and 6.0), and 6.0 will be removed from dict_predictions
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
        self.assertEqual(len(dict_predictions["id1"]), 3)

    def test_load_reference_data_grb_actual_adp(self):
        thematic_dict = {
            "theme_id_1": from_wkt(
                "MultiPolygon (((174184.09476602054201066 171899.68933439542888664, "
                "174400.56834639035514556 171832.959863749332726, "
                "174388.65236948925303295 171770.99678386366576888, "
                "174182.10876987033407204 171836.13745758961886168, "
                "174184.88916448061354458 171873.07698598300339654, "
                "174184.09476602054201066 171899.68933439542888664)))"
            )
        }
        self.sample_aligner.load_thematic_data(DictLoader(thematic_dict))
        # LOAD REFERENCE DICTIONARY
        self.sample_aligner.load_reference_data(
            GRBActualLoader(
                aligner=self.sample_aligner, grb_type=GRBType.ADP, partition=1000
            )
        )
        self.assertGreater(len(self.sample_aligner.dict_reference), 0)

    def test_load_reference_data_grb_actual_gbg(self):
        thematic_dict = {
            "theme_id_1": from_wkt(
                "MultiPolygon (((174184.09476602054201066 171899.68933439542888664, "
                "174400.56834639035514556 171832.959863749332726, "
                "174388.65236948925303295 171770.99678386366576888, "
                "174182.10876987033407204 171836.13745758961886168, "
                "174184.88916448061354458 171873.07698598300339654, "
                "174184.09476602054201066 171899.68933439542888664)))"
            )
        }
        self.sample_aligner.load_thematic_data(DictLoader(thematic_dict))
        # LOAD REFERENCE DICTIONARY
        self.sample_aligner.load_reference_data(
            GRBActualLoader(
                aligner=self.sample_aligner, grb_type=GRBType.GBG, partition=1000
            )
        )
        self.assertGreater(len(self.sample_aligner.dict_reference), 0)

    def test_load_reference_data_grb_actual_knw(self):
        thematic_dict = {
            "theme_id_1": from_wkt(
                "MultiPolygon (((174184.09476602054201066 171899.68933439542888664, "
                "174400.56834639035514556 171832.959863749332726, "
                "174388.65236948925303295 171770.99678386366576888, "
                "174182.10876987033407204 171836.13745758961886168, "
                "174184.88916448061354458 171873.07698598300339654, "
                "174184.09476602054201066 171899.68933439542888664)))"
            )
        }
        self.sample_aligner.load_thematic_data(DictLoader(thematic_dict))
        # LOAD REFERENCE DICTIONARY
        self.sample_aligner.load_reference_data(
            GRBActualLoader(
                aligner=self.sample_aligner, grb_type=GRBType.KNW, partition=1000
            )
        )
        self.sample_aligner.process()
        self.assertGreaterEqual(len(self.sample_aligner.dict_reference), 0)

    def test_all_od_strategies(self):
        # Load thematic data & reference data
        thematic_dict = {
            "theme_id_1": from_wkt("POLYGON ((0 0, 0 9, 5 10, 10 0, 0 0))")
        }
        # ADD A REFERENCE POLYGON TO REFERENCE DICTIONARY
        reference_dict = {"ref_id_1": from_wkt("POLYGON ((0 1, 0 10,8 10,10 1,0 1))")}
        # LOAD THEMATIC DICTIONARY
        self.sample_aligner.load_thematic_data(DictLoader(thematic_dict))
        # LOAD REFERENCE DICTIONARY
        self.sample_aligner.load_reference_data(DictLoader(reference_dict))
        relevant_distance = 1
        for od_strategy in OpenbaarDomeinStrategy:
            process_result = self.sample_aligner.process(
                relevant_distance=relevant_distance,
                od_strategy=od_strategy,
                threshold_overlap_percentage=50,
            )
            self.assertEqual(len(process_result["theme_id_1"][relevant_distance]), 6)

    def test_process_interior_ring(self):
        thematic_dict = {
            "theme_id_1": from_wkt(
                "MultiPolygon (((174370.67910432978533208 171012.2469546866195742, "
                "174461.24052877808571793 171002.71417316573206335, "
                "174429.46459037516615354 170869.25523187351063825, "
                "174373.85669817010057159 170894.6759825958579313, "
                "174369.09030740967136808 170896.6619787460367661, "
                "174370.67910432978533208 171012.2469546866195742),"
                "(174400.07184735251939856 170963.78864862219779752, "
                "174401.26344504262669943 170926.4519209987774957, "
                "174419.53460962430108339 170922.47992869839072227, "
                "174429.06739114515949041 170958.62505863170372322, "
                "174400.07184735251939856 170963.78864862219779752)))"
            )
        }
        self.sample_aligner.load_thematic_data(DictLoader(thematic_dict))
        # LOAD REFERENCE DICTIONARY
        self.sample_aligner.load_reference_data(
            GRBActualLoader(
                aligner=self.sample_aligner, grb_type=GRBType.GBG, partition=1000
            )
        )
        result_dict = self.sample_aligner.process()
        self.assertEqual(len(result_dict), len(thematic_dict))

    def test_process_circle(self):
        # TODO
        geometry = Point(0, 0).buffer(3)
        thematic_dict = {"key": geometry}
        self.sample_aligner.load_thematic_data(DictLoader(thematic_dict))
        # LOAD REFERENCE DICTIONARY
        self.sample_aligner.load_reference_data(
            GRBActualLoader(
                aligner=self.sample_aligner, grb_type=GRBType.GBG, partition=1000
            )
        )
        relevant_distance = 1
        results_dict = self.sample_aligner.process(relevant_distance=relevant_distance)
        self.assertEqual(geometry, results_dict["key"][relevant_distance]["result"])

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
        aligner.dict_thematic, properties, source = thematic_loader.load_data()
        assert aligner.dict_thematic == {"4": shape(geojson["features"][0]["geometry"])}
        self.assertGreater(len(aligner.dict_thematic), 0)

    def test_get_reference_as_geojson(self):
        self.sample_aligner.load_thematic_data(
            DictLoader(
                {"theme_id_1": from_wkt("POLYGON ((0 0, 0 9, 5 10, 10 0, 0 0))")}
            )
        )
        self.sample_aligner.load_reference_data(
            DictLoader({"ref_id_1": from_wkt("POLYGON ((0 1, 0 10,8 10,10 1,0 1))")})
        )
        self.sample_aligner.process()
        self.sample_aligner.get_input_as_geojson()

    def test_fully_aligned_input(self):
        aligned_shape = from_wkt("POLYGON ((0 0, 0 9, 5 10, 10 0, 0 0))")
        loader = DictLoader({"theme_id_1": aligned_shape})
        self.sample_aligner.load_thematic_data(
            DictLoader({"theme_id_1": aligned_shape})
        )
        self.sample_aligner.load_reference_data(DictLoader({"ref_id_1": aligned_shape}))
        relevant_distance = 1
        result = self.sample_aligner.process(relevant_distance=relevant_distance)
        assert result["theme_id_1"][relevant_distance].get("result") == aligned_shape
        assert result["theme_id_1"][relevant_distance].get("result_diff") == Polygon()
        assert (
            result["theme_id_1"][relevant_distance].get("result_diff_min") == Polygon()
        )
        assert (
            result["theme_id_1"][relevant_distance].get("result_diff_plus") == Polygon()
        )

    def test_evaluate(self):
        thematic_dict = {
            "theme_id_1": from_wkt(
                "MultiPolygon (((174180.20077791667426936 171966.14649116666987538, "
                "174415.60530965600628406 171940.9636807945498731, "
                "174388.65236948925303295 171770.99678386366576888, "
                "174182.10876987033407204 171836.13745758961886168, "
                "174184.88916448061354458 171873.07698598300339654, "
                "174180.20077791667426936 171966.14649116666987538)))"
            )
        }
        base_aligner = Aligner()
        base_aligner.load_thematic_data(DictLoader(thematic_dict))
        base_aligner.load_reference_data(
            GRBFiscalParcelLoader(aligner=base_aligner, year="2022", partition=1000)
        )
        relevant_distance = 1
        base_process_result = base_aligner.process(relevant_distance=relevant_distance)
        thematic_dict_formula = {}
        thematic_dict_result = {}
        for key in base_process_result:
            thematic_dict_result[key] = base_process_result[key][relevant_distance][
                "result"
            ]
            thematic_dict_formula[key] = {
                FORMULA_FIELD_NAME: base_aligner.get_brdr_formula(
                    thematic_dict_result[key]
                )
            }
        aligner_result = Aligner()
        aligner_result.load_thematic_data(DictLoader(thematic_dict_result))
        dict_affected, dict_unchanged = get_geoms_affected_by_grb_change(
            aligner=aligner_result,
            grb_type=GRBType.ADP,
            date_start=date(2022, 1, 1),
            date_end=date.today(),
            one_by_one=False,
        )

        actual_aligner = Aligner()
        actual_aligner.load_thematic_data(
            DictLoader(
                data_dict=dict_affected, data_dict_properties=thematic_dict_formula
            )
        )
        loader = GRBActualLoader(
            grb_type=GRBType.ADP, partition=1000, aligner=actual_aligner
        )
        actual_aligner.load_reference_data(loader)
        series = np.arange(0, 200, 10, dtype=int) / 100

        dict_evaluated, prop_dictionary = actual_aligner.compare(
            threshold_area=5,
            threshold_percentage=1,
        )
        fc = get_series_geojson_dict(
            dict_evaluated,
            crs=actual_aligner.CRS,
            id_field=actual_aligner.name_thematic_id,
            series_prop_dict=prop_dictionary,
        )

    def test_fully_aligned_geojson_output(self):
        aligned_shape = from_wkt(
            "MultiPolygon (((173463.11530961000244133 174423.83310307000647299, "
            "173460.22633100001257844 174422.02316300000529736, "
            "173455.24681099998997524 174429.98009100000490434, "
            "173454.4299790000077337 174429.34482699999352917, "
            "173452.06690700000035577 174432.43058700000983663, "
            "173451.25743500000680797 174431.8672589999914635, "
            "173448.74844299998949282 174434.96249100001296028, "
            "173448.5809550000121817 174435.80485899999621324, "
            "173455.39772300000186078 174441.47852299999794923, "
            "173461.44169100001454353 174446.50898700000834651, "
            "173472.15932299999985844 174429.49919500001124106, "
            "173466.18524341000011191 174425.75641125999391079, "
            "173466.18524347001221031 174425.75641117000486702, "
            "173463.11530969000887126 174423.83310300001176074, "
            "173463.11530961000244133 174423.83310307000647299)))"
        )

        self.sample_aligner.load_thematic_data(
            DictLoader({"theme_id_1": aligned_shape})
        )
        self.sample_aligner.load_reference_data(DictLoader({"ref_id_1": aligned_shape}))
        self.sample_aligner.process()
        fcs = self.sample_aligner.get_results_as_geojson(formula=True)
        assert fcs["result"]["features"][0]["properties"]["area"] > 0
        assert fcs["result_diff"]["features"][0]["properties"]["area"] == 0
        assert fcs["result_diff_min"]["features"][0]["properties"]["area"] == 0
        assert fcs["result_diff_plus"]["features"][0]["properties"]["area"] == 0
