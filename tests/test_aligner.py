import json
import os
import unittest
from datetime import date

import numpy as np
from shapely import Point
from shapely import from_wkt
from shapely.geometry import Polygon
from shapely.geometry import shape

from brdr.aligner import Aligner
from brdr.constants import FORMULA_FIELD_NAME, AREA_ATTRIBUTE
from brdr.enums import GRBType, AlignerResultType, Evaluation
from brdr.enums import OpenbaarDomeinStrategy
from brdr.geometry_utils import _grid_bounds, safe_equals
from brdr.geometry_utils import buffer_neg_pos
from brdr.grb import (
    GRBActualLoader,
    GRBFiscalParcelLoader,
    get_affected_by_grb_change,
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
        series = np.arange(0, 310, 10, dtype=int) / 100
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

        series = np.arange(0, 810, 10, dtype=int) / 100
        # predict which relevant distances are interesting to propose as resulting geometry
        dict_series, dict_predictions, diffs = aligner.predictor(
            relevant_distances=series, od_strategy=4, threshold_overlap_percentage=50
        )
        self.assertEqual(len(dict_predictions["id1"]), 3)

    def test_predictor_no_prediction(self):
        """
        Test if no prediction is returned when there are no stable geometries (=no zerostreaks) in the range of relevent_distances
        This testdata has no prediction between 0 and 1 meter, constantly increasing diff (first stable geometry at relevant_distances>2)
        """
        # Initiate an Aligner
        aligner = Aligner()
        # Load thematic data & reference data
        wkt = "Polygon ((174125.66583829143201001 179320.36133541050367057, 174125.66583829294540919 179320.36133543887990527, 174123.56114244047785178 179320.47114562886417843, 174123.60579274909105152 179320.9312379399780184, 174123.65622600910137407 179321.4509197928418871, 174130.00431616476271302 179386.86384879014804028, 174131.20465914410306141 179386.55606853921199217, 174131.47325675669708289 179386.48719735653139651, 174131.4777007331722416 179386.48605787538690493, 174131.48166125148418359 179386.48504236005828716, 174131.56964346842141822 179386.46248281505540945, 174145.93543933291221038 179382.77894541397108696, 174153.73543597472598776 179380.77894627506611869, 174155.22223844396648929 179380.39771487266989425, 174155.31905292189912871 179380.37289064755896106, 174158.08159073328715749 179379.66454761903150938, 174162.78046076380996965 179378.45970914987265132, 174162.77348954943590797 179378.41265345277497545, 174162.76924267943832092 179378.38398708027671091, 174162.76793666195590049 179378.37517146224854514, 174158.29825295542832464 179348.20480644324561581, 174158.06896705977851525 179348.24847994712763466, 174157.82196880917763337 179348.29552723292727023, 174152.4153556079545524 179349.32535831883433275, 174151.38357602406176738 179340.0393457512545865, 174150.45438929076772183 179331.6766651498619467, 174155.32298101272317581 179329.59012584027368575, 174155.97540604253299534 179329.31051516399020329, 174156.06587580699124373 179329.27174254972487688, 174156.3372870393213816 179329.15542325720889494, 174155.77701106332824565 179318.79031770044821315, 174154.96588003114447929 179318.83263741704286076, 174154.95658402171102352 179318.83312241756357253, 174151.97826524809352122 179318.98851313433260657, 174146.5353382924804464 179319.27249193197349086, 174125.66583829143201001 179320.36133541050367057))"

        loader = DictLoader({"id1": from_wkt(wkt)})
        aligner.load_thematic_data(loader)
        loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
        aligner.load_reference_data(loader)

        series = np.arange(0, 110, 10, dtype=int) / 100
        # predict which relevant distances are interesting to propose as resulting geometry
        dict_series, dict_predictions, diffs = aligner.predictor(
            relevant_distances=series, od_strategy=4, threshold_overlap_percentage=50
        )
        self.assertEqual(len(dict_predictions["id1"]), 0)

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
            self.assertEqual(len(process_result["theme_id_1"][relevant_distance]), 7)

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
        geometry = Point(0, 0).buffer(3)
        # geometry = MultiPolygon([geometry])
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
        self.sample_aligner.load_thematic_data(
            DictLoader({"theme_id_1": aligned_shape})
        )
        self.sample_aligner.load_reference_data(DictLoader({"ref_id_1": aligned_shape}))
        relevant_distance = 1
        result = self.sample_aligner.process(relevant_distance=relevant_distance)
        assert safe_equals(
            result["theme_id_1"][relevant_distance].get("result"), aligned_shape
        )
        assert result["theme_id_1"][relevant_distance].get("result_diff").is_empty
        assert result["theme_id_1"][relevant_distance].get("result_diff_min").is_empty
        assert result["theme_id_1"][relevant_distance].get("result_diff_plus").is_empty

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
                FORMULA_FIELD_NAME: json.dumps(
                    base_aligner.get_brdr_formula(thematic_dict_result[key])
                )
            }
            print(key + ": " + thematic_dict_result[key].wkt)
            print(key + ": " + str(thematic_dict_formula[key]))
        base_aligner_result = Aligner()
        base_aligner_result.load_thematic_data(DictLoader(thematic_dict_result))
        affected, unaffected = get_affected_by_grb_change(
            dict_thematic=thematic_dict_result,
            grb_type=GRBType.ADP,
            date_start=date(2022, 1, 1),
            date_end=date.today(),
            one_by_one=False,
            border_distance=relevant_distance,
        )
        if len(affected) == 0:
            print("No affected dicts")
            exit()
        print("Affected_IDs: " + str(affected))
        actual_aligner = Aligner()
        actual_aligner.load_thematic_data(
            DictLoader(
                data_dict=thematic_dict_result,
                data_dict_properties=thematic_dict_formula,
            )
        )
        actual_aligner.load_reference_data(
            GRBActualLoader(
                grb_type=GRBType.ADP, partition=1000, aligner=actual_aligner
            )
        )
        actual_aligner.relevant_distances = np.arange(0, 210, 10, dtype=int) / 100
        dict_evaluated, prop_dictionary = actual_aligner.evaluate(
            ids_to_evaluate=affected, base_formula_field=FORMULA_FIELD_NAME
        )

        fc = get_series_geojson_dict(
            dict_evaluated,
            crs=actual_aligner.CRS,
            id_field=actual_aligner.name_thematic_id,
            series_prop_dict=prop_dictionary,
        )
        print(fc["result"])
        fcs = actual_aligner.get_results_as_geojson(formula=True)

    def test_evaluate_full_parcel_false(self):
        thematic_dict = {
            "theme_id_1": from_wkt(
                "Polygon ((174015.08694170592934825 179025.39916784031083807, 174040.71934808720834553 179031.93985084796440788, 174037.1838437587430235 178986.50862022733781487, 174030.64316075111855753 178982.97311589887249283, 174016.85469387014745735 179002.24161448894301429, 174021.62762471355381422 179005.77711881740833633, 174018.62244603436556645 179008.60552228015149012, 174015.08694170592934825 179025.39916784031083807))"
            )
        }
        aligner = Aligner()
        aligner.load_thematic_data(DictLoader(thematic_dict))
        aligner.load_reference_data(
            GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
        )

        dict_evaluated, prop_dictionary = aligner.evaluate(
            relevant_distances=np.arange(0, 410, 10, dtype=int) / 100, prefer_full=False
        )
        assert len(dict_evaluated["theme_id_1"]) == 1
        assert (
            prop_dictionary["theme_id_1"][3.5]["brdr_evaluation"]
            == Evaluation.TO_CHECK_PREDICTION_MULTI
        )

    def test_evaluate_full_parcel_true(self):
        thematic_dict = {
            "theme_id_1": from_wkt(
                "Polygon ((174015.08694170592934825 179025.39916784031083807, 174040.71934808720834553 179031.93985084796440788, 174037.1838437587430235 178986.50862022733781487, 174030.64316075111855753 178982.97311589887249283, 174016.85469387014745735 179002.24161448894301429, 174021.62762471355381422 179005.77711881740833633, 174018.62244603436556645 179008.60552228015149012, 174015.08694170592934825 179025.39916784031083807))"
            )
        }
        aligner = Aligner()
        aligner.load_thematic_data(DictLoader(thematic_dict))
        aligner.load_reference_data(
            GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
        )

        dict_evaluated, prop_dictionary = aligner.evaluate(
            relevant_distances=np.arange(0, 410, 10, dtype=int) / 100, prefer_full=True
        )
        assert len(dict_evaluated["theme_id_1"]) == 1
        assert (
            prop_dictionary["theme_id_1"][3.5]["brdr_evaluation"]
            == Evaluation.PREDICTION_FULL
        )

    def test_evaluate_all_predictions(self):
        thematic_dict = {
            "theme_id_1": from_wkt(
                "Polygon ((174015.08694170592934825 179025.39916784031083807, 174040.71934808720834553 179031.93985084796440788, 174037.1838437587430235 178986.50862022733781487, 174030.64316075111855753 178982.97311589887249283, 174016.85469387014745735 179002.24161448894301429, 174021.62762471355381422 179005.77711881740833633, 174018.62244603436556645 179008.60552228015149012, 174015.08694170592934825 179025.39916784031083807))"
            )
        }
        aligner = Aligner()
        aligner.load_thematic_data(DictLoader(thematic_dict))
        aligner.load_reference_data(
            GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
        )

        dict_evaluated, prop_dictionary = aligner.evaluate(
            relevant_distances=np.arange(0, 410, 10, dtype=int) / 100,
            prefer_full=True,
            all_predictions=True,
        )
        assert len(dict_evaluated["theme_id_1"]) == 2
        assert (
            prop_dictionary["theme_id_1"][3.5]["brdr_evaluation"]
            == Evaluation.PREDICTION_FULL
        )

    def test_remark_for_poly_multipoly(self):
        shape = from_wkt(
            "MultiPolygon(((48893.03662109375 214362.93756103515625, 48890.8258056640625 214368.482666015625, 48890.7159423828125 214368.44110107421875, 48887.6488037109375 214367.2845458984375, 48886.3800048828125 214368.68017578125, 48885.1068115234375 214370.08062744140625, 48884.3330078125 214369.782470703125, 48882.563720703125 214369.10064697265625, 48882.1116943359375 214370.1346435546875, 48878.5626220703125 214368.70196533203125, 48877.839111328125 214368.40997314453125, 48877.2352294921875 214369.79376220703125, 48876.7911376953125 214369.60687255859375, 48875.0850830078125 214373.62353515625, 48875.478759765625 214373.8182373046875, 48881.5286865234375 214376.81109619140625, 48885.10546875 214372.36151123046875, 48887.0050048828125 214370.08538818359375, 48888.4698486328125 214368.330078125, 48890.366943359375 214369.2685546875, 48901.0638427734375 214374.56024169921875, 48905.0159912109375 214369.61175537109375, 48904.472900390625 214367.53851318359375, 48893.03662109375 214362.93756103515625)))"
        )
        self.sample_aligner.load_thematic_data(DictLoader({"theme_id_1": shape}))
        self.sample_aligner.load_reference_data(
            GRBActualLoader(
                grb_type=GRBType.ADP, partition=1000, aligner=self.sample_aligner
            )
        )
        self.sample_aligner.process(relevant_distances=[2])
        assert self.sample_aligner.dict_processresults["theme_id_1"][2]["remark"] != ""

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
        assert fcs["result"]["features"][0]["properties"][AREA_ATTRIBUTE] > 0
        assert fcs["result_diff"]["features"][0]["properties"][AREA_ATTRIBUTE] == 0
        assert fcs["result_diff_min"]["features"][0]["properties"][AREA_ATTRIBUTE] == 0
        assert fcs["result_diff_plus"]["features"][0]["properties"][AREA_ATTRIBUTE] == 0
