import unittest
from types import SimpleNamespace
from unittest.mock import patch

from shapely.geometry import Point

from brdr.aligner import AlignerResult
from brdr.constants import (
    EVALUATION_FIELD_NAME,
    PREDICTION_COUNT,
    PREDICTION_SCORE,
    STABILITY,
)
from brdr.enums import Evaluation
from brdr.evaluator import AlignerEvaluator, ConservativeAlignerEvaluator
from brdr.predictor import AlignerPredictor, ConservativeAlignerPredictor


class TestConservativePredictor(unittest.TestCase):
    def test_keeps_only_high_confidence_stable_predictions(self):
        process_results = {
            "theme_1": {
                0.0: {
                    "result": Point(0, 0),
                    "properties": {PREDICTION_SCORE: 42, STABILITY: True},
                },
                0.5: {
                    "result": Point(0, 0),
                    "properties": {PREDICTION_SCORE: 80, STABILITY: True},
                },
                1.0: {
                    "result": Point(0, 0),
                    "properties": {PREDICTION_SCORE: 78, STABILITY: True},
                },
                1.5: {
                    "result": Point(0, 0),
                    "properties": {PREDICTION_SCORE: 90, STABILITY: False},
                },
            }
        }
        base_result = AlignerResult(process_results)
        predictor = ConservativeAlignerPredictor(
            min_prediction_score=70,
            score_tolerance_from_best=1,
            require_stable=True,
        )

        with patch.object(AlignerPredictor, "predict", return_value=base_result):
            out = predictor.predict(
                aligner=SimpleNamespace(), relevant_distances=[0.0, 0.5, 1.0, 1.5]
            )

        props_05 = out.results["theme_1"][0.5]["properties"]
        props_10 = out.results["theme_1"][1.0]["properties"]
        props_00 = out.results["theme_1"][0.0]["properties"]
        props_15 = out.results["theme_1"][1.5]["properties"]

        self.assertIn(PREDICTION_SCORE, props_05)
        self.assertNotIn(PREDICTION_SCORE, props_10)
        self.assertNotIn(PREDICTION_SCORE, props_00)
        self.assertNotIn(PREDICTION_SCORE, props_15)
        self.assertEqual(props_05[PREDICTION_COUNT], 1)


class TestConservativeEvaluator(unittest.TestCase):
    def test_force_single_prediction_keeps_best(self):
        evaluated = {
            "theme_1": {
                0.0: {
                    "result": Point(0, 0),
                    "properties": {
                        PREDICTION_SCORE: 71,
                        EVALUATION_FIELD_NAME: Evaluation.TO_CHECK_PREDICTION_MULTI,
                    },
                },
                0.5: {
                    "result": Point(0, 0),
                    "properties": {
                        PREDICTION_SCORE: 90,
                        EVALUATION_FIELD_NAME: Evaluation.TO_CHECK_PREDICTION_MULTI,
                    },
                },
                1.0: {
                    "result": Point(0, 0),
                    "properties": {
                        PREDICTION_SCORE: 80,
                        EVALUATION_FIELD_NAME: Evaluation.TO_CHECK_PREDICTION_MULTI,
                    },
                },
            }
        }
        evaluator = ConservativeAlignerEvaluator(
            ambiguity_delta=5, force_single_prediction=True
        )
        aligner = SimpleNamespace(
            thematic_data=SimpleNamespace(
                features={"theme_1": SimpleNamespace(geometry=Point(0, 0))}
            )
        )

        with patch.object(
            AlignerEvaluator, "evaluate", return_value=AlignerResult(evaluated)
        ):
            out = evaluator.evaluate(
                aligner=aligner, relevant_distances=[0.0, 0.5, 1.0]
            )

        self.assertEqual(list(out.results["theme_1"].keys()), [0.5])
        self.assertEqual(out.results["theme_1"][0.5]["properties"][PREDICTION_COUNT], 1)

    def test_ambiguity_returns_original(self):
        evaluated = {
            "theme_1": {
                0.0: {
                    "result": Point(0, 0),
                    "properties": {
                        PREDICTION_SCORE: 70,
                        EVALUATION_FIELD_NAME: Evaluation.TO_CHECK_PREDICTION_MULTI,
                    },
                },
                0.5: {
                    "result": Point(0, 0),
                    "properties": {
                        PREDICTION_SCORE: 90,
                        EVALUATION_FIELD_NAME: Evaluation.TO_CHECK_PREDICTION_MULTI,
                    },
                },
                1.0: {
                    "result": Point(0, 0),
                    "properties": {
                        PREDICTION_SCORE: 88,
                        EVALUATION_FIELD_NAME: Evaluation.TO_CHECK_PREDICTION_MULTI,
                    },
                },
            }
        }
        evaluator = ConservativeAlignerEvaluator(
            ambiguity_delta=5, force_single_prediction=True
        )
        aligner = SimpleNamespace(
            thematic_data=SimpleNamespace(
                features={"theme_1": SimpleNamespace(geometry=Point(1, 1))}
            )
        )

        def fake_update(**kwargs):
            process_results = kwargs["process_results_evaluated"]
            process_results["theme_1"][0] = {
                "result": Point(1, 1),
                "properties": {
                    EVALUATION_FIELD_NAME: Evaluation.TO_CHECK_ORIGINAL,
                    PREDICTION_SCORE: -1,
                },
            }
            return process_results

        with patch.object(
            AlignerEvaluator, "evaluate", return_value=AlignerResult(evaluated)
        ):
            with patch.object(
                ConservativeAlignerEvaluator,
                "update_evaluation_with_original",
                side_effect=fake_update,
            ):
                out = evaluator.evaluate(
                    aligner=aligner, relevant_distances=[0.0, 0.5, 1.0]
                )

        self.assertEqual(list(out.results["theme_1"].keys()), [0])
        self.assertEqual(
            out.results["theme_1"][0]["properties"][EVALUATION_FIELD_NAME],
            Evaluation.TO_CHECK_ORIGINAL,
        )
