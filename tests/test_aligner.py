import os
import unittest

import numpy as np
from shapely import Point
from shapely import from_wkt
from shapely.geometry import Polygon
from shapely.geometry import shape

from brdr.aligner import Aligner
from brdr.constants import AREA_ATTRIBUTE
from brdr.enums import GRBType, AlignerResultType
from brdr.enums import OpenDomainStrategy
from brdr.geometry_utils import _grid_bounds, safe_equals, geom_from_wkt
from brdr.geometry_utils import buffer_neg_pos
from brdr.grb import (
    GRBActualLoader,
)
from brdr.loader import GeoJsonLoader, DictLoader
from brdr.typings import FeatureCollection, ProcessResult


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
        self.assertEqual(len(dict_predictions["id1"]), 2)

    def test_line(self):
        aligner = Aligner(max_workers=-1)
        wkt = "MULTILINESTRING ((174024.1298775521281641 179420.42107488788315095, 174042.22504128722357564 179413.1830093938333448, 174044.36356063772109337 179418.28255553735652938, 174049.29860529277357273 179418.94056149138486944, 174050.28561422377242707 179422.72409572688047774, 174054.39815143629675731 179421.90158828438143246, 174057.52367971779312938 179420.75007786488276906, 174054.56265292479656637 179408.08346325031016022, 174047.16008594224695116 179398.54237691726302728, 174036.960993655200582 179404.29992901478544809, 174032.84845644267625175 179400.51639477928983979, 174029.88742964965058491 179397.88437096326379105))"
        thematic_dict = {"theme_id_1": from_wkt(wkt)}

        loader = DictLoader(thematic_dict)
        aligner.load_thematic_data(loader)
        loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
        aligner.load_reference_data(loader)
        relevant_distance = 3
        dict_processresults = aligner.process(relevant_distance=relevant_distance)
        self.assertEqual(
            dict_processresults["theme_id_1"][relevant_distance]["result"].geom_type,
            "LineString",
        )

    def test_reference_mix(self):
        "reference exists out of points, lines and polygons"
        aligner = Aligner(max_workers=-1)
        wkt = "POLYGON ((174043.1035556931165047 179299.26716663906699978, 174042.53709723605425097 179306.22651339753065258, 174048.52537235378986225 179305.9837454873486422, 174055.16102856534416787 179305.49820966698462144, 174055.16102856534416787 179293.19796888460405171, 174071.26463327385135926 179293.92627261512097903, 174071.34555591057869606 179287.69522958720335737, 174083.40302878280635923 179282.83987138362135738, 174072.154782277852064 179278.14635845352313481, 174072.72124073494342156 179260.58614628392388113, 174056.69855866313446313 179258.88677091267891228, 174048.60629499051719904 179273.21007761321379803, 174038.97650122008053586 179273.85745870703249238, 174041.97063877896289341 179286.72415794650441967, 174044.96477633781614713 179292.71243306424003094, 174035.25405993068125099 179292.14597460714867339, 174043.1035556931165047 179299.26716663906699978))"
        thematic_dict = {"theme_id_1": from_wkt(wkt)}
        loader = DictLoader(thematic_dict)
        aligner.load_thematic_data(loader)
        point_1 = geom_from_wkt(
            "POINT (174043.75093678693519905 179298.98393741055042483)"
        )
        point_2 = geom_from_wkt(
            "POINT (174079.35689694649772719 179283.04217797549790703)"
        )
        line_1 = geom_from_wkt(
            "LINESTRING (174021.41628905048128217 179291.53905483175185509, 174081.29904022792470641 179294.45226975387777202)"
        )
        line_2 = geom_from_wkt(
            "LINESTRING (174039.13834649353520945 179282.79941006531589665, 174038.16727485283627175 179272.92684838469722308, 174052.65242682682583109 179273.33146156833390705, 174093.43743573685060255 179285.46985707725980319)"
        )
        polygon_1 = geom_from_wkt(
            "POLYGON ((174029.67039799655321985 179306.67158789956010878, 174048.44444971703342162 179306.18605207919608802, 174048.92998553739744239 179258.92723223107168451, 174030.15593381691724062 179259.41276805143570527, 174029.67039799655321985 179306.67158789956010878))"
        )
        polygon_2 = geom_from_wkt(
            "POLYGON ((174054.91826065513305366 179305.53867098540649749, 174056.37486811622511595 179258.76538695761701092, 174072.72124073491431773 179260.38383969213464297, 174070.61725218003266491 179306.99527844646945596, 174054.91826065513305366 179305.53867098540649749))"
        )
        reference_dict = {
            "ref_id_1": point_1,
            "ref_id_2": point_2,
            "ref_id_3": line_1,
            "ref_id_4": line_2,
            "ref_id_5": polygon_1,
            "ref_id_6": polygon_2,
        }

        loader = DictLoader(reference_dict)
        aligner.load_reference_data(loader)
        relevant_distance = 5
        dict_processresults = aligner.process(relevant_distance=relevant_distance)
        self.assertEqual(
            dict_processresults["theme_id_1"][relevant_distance]["result"].geom_type,
            "Polygon",
        )

    def test_predictor_no_prediction(self):
        """
        Test if no prediction is returned when there are no stable geometries (=no zerostreaks) in the range of relevent_distances
        This testdata has no prediction between 0 and 1 meter, constantly increasing diff (first stable geometry at relevant_distances>2)
        """
        # Initiate an Aligner
        aligner = Aligner(max_workers=-1)
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
        for od_strategy in OpenDomainStrategy:
            process_result = self.sample_aligner.process(
                relevant_distance=relevant_distance,
                od_strategy=od_strategy,
                threshold_overlap_percentage=50,
            )
            self.assertEqual(len(process_result["theme_id_1"][relevant_distance]), 8)

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

    def test_remark_for_poly_multipoly(self):
        # TODO correct test
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
        assert self.sample_aligner.dict_processresults["theme_id_1"][2]["remark"] == ""

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
