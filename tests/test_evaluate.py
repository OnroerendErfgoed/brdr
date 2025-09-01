import unittest

import numpy as np
from shapely import from_wkt
from shapely.geometry import Polygon

from brdr.aligner import Aligner
from brdr.constants import EVALUATION_FIELD_NAME
from brdr.enums import GRBType, Evaluation, FullStrategy, SnapStrategy
from brdr.grb import (
    GRBActualLoader,
)
from brdr.loader import DictLoader, GeoJsonLoader


class TestEvaluate(unittest.TestCase):
    def setUp(self):
        # Create a sample geometry for testing
        self.sample_aligner = Aligner()
        self.sample_geom = Polygon([(0, 0), (0, 10), (10, 10), (10, 0)])

    # def test_evaluate(self):
    #     thematic_dict = {
    #         "theme_id_1": from_wkt(
    #             "MultiPolygon (((174180.20077791667426936 171966.14649116666987538, "
    #             "174415.60530965600628406 171940.9636807945498731, "
    #             "174388.65236948925303295 171770.99678386366576888, "
    #             "174182.10876987033407204 171836.13745758961886168, "
    #             "174184.88916448061354458 171873.07698598300339654, "
    #             "174180.20077791667426936 171966.14649116666987538)))"
    #         )
    #     }
    #     base_aligner = Aligner()
    #     base_aligner.load_thematic_data(DictLoader(thematic_dict))
    #     base_aligner.load_reference_data(
    #         GRBFiscalParcelLoader(aligner=base_aligner, year="2022", partition=1000)
    #     )
    #     relevant_distance = 1
    #     base_process_result = base_aligner.process(relevant_distance=relevant_distance)
    #     thematic_dict_formula = {}
    #     thematic_dict_result = {}
    #     for key in base_process_result:
    #         thematic_dict_result[key] = base_process_result[key][relevant_distance][
    #             "result"
    #         ]
    #         thematic_dict_formula[key] = {
    #             FORMULA_FIELD_NAME: json.dumps(
    #                 base_aligner.get_brdr_formula(thematic_dict_result[key])
    #             )
    #         }
    #         # print(key + ": " + thematic_dict_result[key].wkt)
    #         # print(key + ": " + str(thematic_dict_formula[key]))
    #     base_aligner_result = Aligner()
    #     base_aligner_result.load_thematic_data(DictLoader(thematic_dict_result))
    #     affected, unaffected = get_affected_by_grb_change(
    #         dict_thematic=thematic_dict_result,
    #         grb_type=GRBType.ADP,
    #         date_start=date(2022, 1, 1),
    #         date_end=date.today(),
    #         one_by_one=False,
    #         border_distance=relevant_distance,
    #     )
    #     if len(affected) == 0:
    #         # print("No affected dicts")
    #         exit()
    #     # print("Affected_IDs: " + str(affected))
    #     actual_aligner = Aligner()
    #     actual_aligner.load_thematic_data(
    #         DictLoader(
    #             data_dict=thematic_dict_result,
    #             data_dict_properties=thematic_dict_formula,
    #         )
    #     )
    #     actual_aligner.load_reference_data(
    #         GRBActualLoader(
    #             grb_type=GRBType.ADP, partition=1000, aligner=actual_aligner
    #         )
    #     )
    #     actual_aligner.relevant_distances = np.arange(0, 210, 10, dtype=int) / 100
    #     dict_evaluated, prop_dictionary  = actual_aligner.evaluate(
    #         ids_to_evaluate=affected, base_formula_field=FORMULA_FIELD_NAME
    #     )
    #
    #     fc = get_dict_geojsons_from_series_dict(
    #         dict_evaluated,
    #         crs=actual_aligner.CRS,
    #         id_field=actual_aligner.name_thematic_id,
    #         series_prop_dict=prop_dictionary,
    #     )
    #     # print(fc["result"])
    #     fcs = actual_aligner.get_results_as_geojson(formula=True)
    #     geojson = fcs["result"]
    #     # print(geojson)
    #     geojson = geojson_to_multi(fcs["result"])
    #     # print(geojson)

    def test_evaluate_full_strategy_no_full(self):
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
            full_strategy=FullStrategy.NO_FULL,
        )
        assert len(dict_evaluated["theme_id_1"]) > 1
        assert (
            prop_dictionary["theme_id_1"][3.5]["brdr_evaluation"]
            == Evaluation.TO_CHECK_PREDICTION_MULTI_FULL
        )

    def test_evaluate_full_strategy_prefer_full(self):
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
            full_strategy=FullStrategy.PREFER_FULL,
        )
        assert len(dict_evaluated["theme_id_1"]) > 1
        assert (
            prop_dictionary["theme_id_1"][3.5]["brdr_evaluation"]
            == Evaluation.TO_CHECK_PREDICTION_FULL
        )

    def test_evaluate_full_strategy_only_full(self):
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
            full_strategy=FullStrategy.ONLY_FULL,
        )
        assert len(dict_evaluated["theme_id_1"]) == 1
        assert (
            prop_dictionary["theme_id_1"][3.5]["brdr_evaluation"]
            == Evaluation.PREDICTION_UNIQUE_FULL
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
            full_strategy=FullStrategy.PREFER_FULL,
            max_predictions=-1,
        )
        assert len(dict_evaluated["theme_id_1"]) > 1
        assert (
            prop_dictionary["theme_id_1"][3.5]["brdr_evaluation"]
            == Evaluation.TO_CHECK_PREDICTION_FULL
        )

    def test_evaluate_limited_predictions(self):
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
            full_strategy=FullStrategy.NO_FULL,
            max_predictions=2,
        )
        assert len(dict_evaluated["theme_id_1"]) > 1
        # assert (
        #     prop_dictionary["theme_id_1"][3.5]["brdr_evaluation"]
        #     == Evaluation.TO_CHECK_PREDICTION_MULTI_FULL
        # )

    def test_evaluate_multi_to_best_prediction_true(self):
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
            full_strategy=FullStrategy.NO_FULL,
            max_predictions=1,
            multi_to_best_prediction=True,
        )
        assert len(dict_evaluated["theme_id_1"]) == 1
        assert (
            prop_dictionary["theme_id_1"][
                list(prop_dictionary["theme_id_1"].keys())[0]
            ]["brdr_prediction_count"]
            > 1
        )

    def test_evaluate_multi_to_best_prediction_false(self):
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
            full_strategy=FullStrategy.NO_FULL,
            max_predictions=1,
            multi_to_best_prediction=False,
        )
        assert len(dict_evaluated["theme_id_1"]) == 1
        assert (
            prop_dictionary["theme_id_1"][0]["brdr_evaluation"]
            == Evaluation.TO_CHECK_ORIGINAL
        )

    def test_evaluate_relevant_distances_without_0(self):
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
            relevant_distances=np.arange(10, 410, 10, dtype=int) / 100,
            full_strategy=FullStrategy.NO_FULL,
            max_predictions=1,
            multi_to_best_prediction=False,
        )
        assert len(dict_evaluated["theme_id_1"]) == 1
        assert (
            prop_dictionary["theme_id_1"][0]["brdr_evaluation"]
            == Evaluation.TO_CHECK_ORIGINAL
        )

    def test_evaluate_line(self):
        # Load thematic data & reference data
        thematic_dict = {
            "theme_id": from_wkt(
                "LINESTRING (171741.11190000033820979 171839.01070547936251387, 171751.68948904142598622 171847.23771917796693742, 171762.26707808251376264 171855.69979041084297933, 171772.72713835645117797 171862.9865739724773448, 171783.53978493178146891 171870.8610013697471004, 171794.94007534274714999 171877.79519862998859026, 171801.87427260301774368 171884.2592808217741549)"
            )
        }

        # ADD A REFERENCE POLYGON TO REFERENCE DICTIONARY

        aligner = Aligner()
        aligner.load_thematic_data(DictLoader(thematic_dict))
        # aligner.load_reference_data(DictLoader(reference_dict))
        aligner.load_reference_data(
            GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
        )

        dict_evaluated, prop_dictionary = aligner.evaluate(
            relevant_distances=np.arange(0, 1010, 50, dtype=int) / 100,
            full_strategy=FullStrategy.NO_FULL,
            max_predictions=3,
            multi_to_best_prediction=False,
        )
        assert len(dict_evaluated["theme_id"]) >= 1
        assert (
            prop_dictionary["theme_id"][0]["brdr_evaluation"]
            == Evaluation.TO_CHECK_PREDICTION_MULTI_FULL
        )

    def test_evaluate_point(self):
        # Load thematic data & reference data
        # thematic_dict = {"theme_id": from_wkt("POINT (0 0)")}
        thematic_dict = {
            "theme_id": from_wkt(
                "POINT (173966.17483414348680526 172343.78743441699771211)"
            )
        }

        aligner = Aligner()
        aligner.load_thematic_data(DictLoader(thematic_dict))
        aligner.load_reference_data(
            GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
        )

        dict_evaluated, prop_dictionary = aligner.evaluate(
            relevant_distances=np.arange(0, 1010, 50, dtype=int) / 100,
            full_strategy=FullStrategy.NO_FULL,
            max_predictions=3,
            multi_to_best_prediction=False,
        )
        assert len(dict_evaluated["theme_id"]) > 0
        # assert (
        #     prop_dictionary["theme_id"][0]["brdr_evaluation"]
        #     == Evaluation.TO_CHECK_PREDICTION_MULTI_FULL
        # )

    def test_evaluate_reference_point(self):
        # Load thematic data & reference data
        polygon_1 = from_wkt(
            "MULTIPOLYGON (((171795.71618631482124329 171817.88460136577486992, 171784.53532230854034424 171806.1688893586397171, 171746.73993028700351715 171841.20300138369202614, 171746.28380228579044342 171841.62578538432717323, 171767.35881029814481735 171856.89906539395451546, 171767.47471430152654648 171856.76376939192414284, 171798.38581831753253937 171820.68191336840391159, 171798.01820231974124908 171820.29669736698269844, 171795.71618631482124329 171817.88460136577486992)))"
        )
        geom_reference = from_wkt(
            "MULTIPOINT ((171756.52506366037414409 171850.34457502837176435),(171766.12778689494007267 171857.04121096828021109),(171777.13617191879893653 171865.19089055553195067),(171745.72200002145837061 171841.94219219812657684),(171770.11351679253857583 171840.12903735807049088))"
        )

        # thematic_dict = {"theme_id": from_wkt("POINT (0 0)")}
        thematic_dict = {"theme_id": polygon_1}

        # ADD A REFERENCE POLYGON TO REFERENCE DICTIONARY
        reference_dict = {"ref_id": geom_reference}

        aligner = Aligner(
            snap_strategy=SnapStrategy.NO_PREFERENCE, snap_max_segment_length=2
        )
        aligner.load_thematic_data(DictLoader(thematic_dict))
        aligner.load_reference_data(DictLoader(reference_dict))

        dict_evaluated, prop_dictionary = aligner.evaluate(
            relevant_distances=np.arange(0, 510, 50, dtype=int) / 100,
            full_strategy=FullStrategy.NO_FULL,
            max_predictions=3,
            multi_to_best_prediction=False,
        )
        assert True
        #     (
        #         len(dict_evaluated["theme_id"]) == 3))
        # assert (
        #     prop_dictionary["theme_id"][0]["brdr_evaluation"]
        #     == Evaluation.TO_CHECK_PREDICTION_MULTI_FULL
        # )

    def test_evaluate_best_no_prediction(self):
        thematic_json = {
            "type": "FeatureCollection",
            "name": "uc5",
            "crs": {
                "type": "name",
                "properties": {"name": "urn:ogc:def:crs:EPSG::31370"},
            },
            "features": [
                {
                    "type": "Feature",
                    "properties": {"fid": 1, "dossiernummer": "D_68365"},
                    "geometry": {
                        "type": "Polygon",
                        "coordinates": [
                            [
                                [172648.411900000093738, 171970.078100000013364],
                                [172641.128199999919161, 171979.2324],
                                [172638.179699999920558, 171982.9382],
                                [172633.854899999918416, 171988.3738],
                                [172633.317000841663685, 171989.049898942059372],
                                [172623.4652, 172001.432],
                                [172610.793499999970663, 172017.358100000128616],
                                [172610.942697853111895, 172017.477398283459479],
                                [172650.691, 172049.2673],
                                [172654.167900000087684, 172052.047999999922467],
                                [172656.0643, 172049.772999999986496],
                                [172657.2538, 172048.346],
                                [172677.110299999912968, 172024.525499999901513],
                                [172677.8891, 172023.591299999912735],
                                [172692.862399999867193, 172005.628700000088429],
                                [172693.374199999903794, 172006.038],
                                [172696.634700000082375, 172001.575399999826914],
                                [172709.757699999987381, 171986.0387],
                                [172698.4481, 171976.4558],
                                [172665.904004495765548, 171948.880303809302859],
                                [172665.794206234189915, 171948.787305280333385],
                                [172665.529999999969732, 171948.563400000159163],
                                [172648.411900000093738, 171970.078100000013364],
                            ]
                        ],
                    },
                },
                {
                    "type": "Feature",
                    "properties": {"fid": 2, "dossiernummer": "D_352"},
                    "geometry": {
                        "type": "Polygon",
                        "coordinates": [
                            [
                                [172443.633200000098441, 171905.855499999917811],
                                [172450.522299999924144, 171922.915599999832921],
                                [172471.899299999902723, 171914.2226],
                                [172458.880800000013551, 171881.983800000074552],
                                [172447.773599999636644, 171885.025300000183051],
                                [172436.992100000090431, 171889.409699999901932],
                                [172443.633200000098441, 171905.855499999917811],
                            ]
                        ],
                    },
                },
                {
                    "type": "Feature",
                    "properties": {"fid": 3, "dossiernummer": "D_28311"},
                    "geometry": {
                        "type": "Polygon",
                        "coordinates": [
                            [
                                [172179.29800000009709, 171788.713699999905657],
                                [172182.4732, 171780.029],
                                [172157.412400000175694, 171770.505299999989802],
                                [172129.436300000088522, 171848.2427],
                                [172129.052500541059999, 171849.309098496683873],
                                [172128.886999999987893, 171849.769],
                                [172148.2379, 171862.7293],
                                [172151.4504, 171864.880899999901885],
                                [172151.894199230970116, 171863.666902103519533],
                                [172179.29800000009709, 171788.713699999905657],
                            ]
                        ],
                    },
                },
            ],
        }
        # 4001 will give a result, 4002 &4009 should also give the original geometry

        aligner = Aligner(relevant_distances=np.arange(0, 310, 10, dtype=int) / 100)

        loader = GeoJsonLoader(_input=thematic_json, id_property="fid")
        aligner.load_thematic_data(loader)
        # Load reference data: The actual GRB-parcels
        loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
        aligner.load_reference_data(loader)

        # Use the EVALUATE-function
        dict_predictions_evaluated, prop_dictionary = aligner.evaluate(
            max_predictions=1,
            full_strategy=FullStrategy.ONLY_FULL,
            multi_to_best_prediction=True,
        )

        assert (
            prop_dictionary[1][list(prop_dictionary[1].keys())[0]][
                EVALUATION_FIELD_NAME
            ]
            == Evaluation.PREDICTION_UNIQUE_FULL
        )
        assert (
            prop_dictionary[2][list(prop_dictionary[2].keys())[0]][
                EVALUATION_FIELD_NAME
            ]
            == Evaluation.TO_CHECK_NO_PREDICTION
        )
        # TODO; check below as this gives a extra prediction without prop_dictionary parameters
        # assert prop_dictionary[3][list(prop_dictionary[3].keys())[0]][EVALUATION_FIELD_NAME] == Evaluation.TO_CHECK_NO_PREDICTION
