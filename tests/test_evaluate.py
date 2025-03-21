import json
import unittest
from datetime import date

import numpy as np
from shapely import from_wkt
from shapely.geometry import Polygon

from brdr.aligner import Aligner
from brdr.constants import FORMULA_FIELD_NAME
from brdr.enums import GRBType, Evaluation, FullStrategy
from brdr.grb import (
    GRBActualLoader,
    GRBFiscalParcelLoader,
    get_affected_by_grb_change,
)
from brdr.loader import DictLoader
from brdr.utils import get_series_geojson_dict


class TestEvaluate(unittest.TestCase):
    def setUp(self):
        # Create a sample geometry for testing
        self.sample_aligner = Aligner()
        self.sample_geom = Polygon([(0, 0), (0, 10), (10, 10), (10, 0)])

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
        assert len(dict_evaluated["theme_id_1"]) == 3
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
        assert len(dict_evaluated["theme_id_1"]) == 3
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
        assert len(dict_evaluated["theme_id_1"]) == 3
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
        assert len(dict_evaluated["theme_id_1"]) == 2
        assert (
            prop_dictionary["theme_id_1"][3.5]["brdr_evaluation"]
            == Evaluation.TO_CHECK_PREDICTION_MULTI_FULL
        )

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
            prop_dictionary["theme_id_1"][0.7]["brdr_evaluation"]
            == Evaluation.TO_CHECK_PREDICTION_MULTI
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
