import json
from abc import ABC
from abc import abstractmethod
from copy import deepcopy
from dataclasses import dataclass
from typing import Any
from typing import Iterable
from typing import List
from typing import Optional
from typing import TYPE_CHECKING

import numpy as np
from shapely.geometry.base import BaseGeometry

from brdr.constants import DIFF_AREA_FIELD_NAME
from brdr.constants import DIFF_PERCENTAGE_FIELD_NAME
from brdr.constants import EQUAL_REFERENCE_FEATURES_FIELD_NAME
from brdr.constants import EVALUATION_FIELD_NAME
from brdr.constants import FULL_ACTUAL_FIELD_NAME
from brdr.constants import FULL_BASE_FIELD_NAME
from brdr.constants import METADATA_FIELD_NAME
from brdr.constants import OD_ALIKE_FIELD_NAME
from brdr.constants import PREDICTION_COUNT
from brdr.constants import PREDICTION_SCORE
from brdr.constants import RELEVANT_DISTANCE_DECIMALS
from brdr.constants import REMARK_FIELD_NAME
from brdr.enums import AlignerResultType
from brdr.enums import Evaluation
from brdr.enums import FullReferenceStrategy
from brdr.enums import ProcessRemark
from brdr.logger import Logger
from brdr.metadata import (
    reverse_metadata_observations_to_brdr_observation as _reverse_metadata_observations_core,
)
from brdr.typings import InputId
from brdr.utils import deep_merge
from brdr.utils import is_brdr_observation

if TYPE_CHECKING:
    from brdr.aligner import Aligner
    from brdr.aligner import AlignerResult


def _reverse_metadata_observations_to_brdr_observation(metadata: dict) -> dict:
    return _reverse_metadata_observations_core(metadata)


@dataclass
class _EvaluationRuntime:
    process_results_predictions: dict
    process_results_temp_predictions: dict
    process_results_evaluated: dict
    metadata_field: str
    full_reference_strategy: FullReferenceStrategy
    max_predictions: int
    multi_to_best_prediction: bool


class BaseEvaluator(ABC):
    """
    Abstract base class for evaluation strategies.

    An evaluator compares predictions with the source observation context and
    selects the most relevant candidate geometries.
    """

    def __init__(self, feedback=None):
        self.logger = Logger(feedback)

    @abstractmethod
    def evaluate(
        self,
        *,
        aligner: "Aligner",
        relevant_distances: Optional[Iterable[float]] = None,
        thematic_ids: Optional[List[InputId]] = None,
        metadata_field: str = METADATA_FIELD_NAME,
        full_reference_strategy: FullReferenceStrategy = FullReferenceStrategy.NO_FULL_REFERENCE,
        max_predictions: int = -1,
        multi_to_best_prediction: bool = True,
    ) -> "AlignerResult":
        """Execute evaluation strategy for the given aligner context."""
        pass

    @abstractmethod
    def update_evaluation_with_original(
        self,
        *,
        aligner: "Aligner",
        metadata_field: str,
        original_geometry: BaseGeometry,
        process_results_evaluated,
        theme_id: str | int,
        evaluation: Evaluation,
    ) -> Any:
        """Update evaluated output by returning the original geometry."""
        pass

    @abstractmethod
    def get_observation_properties(
        self,
        *,
        aligner: "Aligner",
        process_result,
    ) -> dict:
        """Build observation properties for a process result."""
        pass

    @abstractmethod
    def get_observation_comparison_properties(
        self,
        *,
        aligner: "Aligner",
        process_result,
        base_brdr_observation=None,
    ) -> dict:
        """Compare current observation against base observation."""
        pass

    @abstractmethod
    def get_brdr_observation_from_properties(
        self,
        *,
        aligner: "Aligner",
        id_theme: Any,
        base_metadata_field: str,
    ) -> dict:
        """Extract a BRDR observation from thematic feature properties."""
        pass


class AlignerEvaluator(BaseEvaluator):
    """
    Default evaluator implementation with the original Aligner evaluate logic.
    """

    def evaluate(
        self,
        *,
        aligner: "Aligner",
        relevant_distances: Optional[Iterable[float]] = None,
        thematic_ids: Optional[List[InputId]] = None,
        metadata_field: str = METADATA_FIELD_NAME,
        full_reference_strategy: FullReferenceStrategy = FullReferenceStrategy.NO_FULL_REFERENCE,
        max_predictions: int = -1,
        multi_to_best_prediction: bool = True,
    ) -> "AlignerResult":
        from brdr.aligner import AlignerResult

        thematic_ids, calculate_zeros = self._resolve_thematic_ids(
            aligner=aligner, thematic_ids=thematic_ids
        )
        relevant_distances = self._resolve_relevant_distances(relevant_distances)
        self._validate_evaluation_inputs(
            aligner=aligner,
            thematic_ids=thematic_ids,
            relevant_distances=relevant_distances,
        )

        process_results, process_results_predictions = self._prepare_process_results(
            aligner=aligner,
            thematic_ids=thematic_ids,
            relevant_distances=relevant_distances,
            calculate_zeros=calculate_zeros,
        )
        runtime = _EvaluationRuntime(
            process_results_predictions=process_results_predictions,
            process_results_temp_predictions=deepcopy(process_results_predictions),
            process_results_evaluated=deepcopy(process_results),
            metadata_field=metadata_field,
            full_reference_strategy=full_reference_strategy,
            max_predictions=max_predictions,
            multi_to_best_prediction=multi_to_best_prediction,
        )

        for theme_id, feat in aligner.thematic_data.features.items():
            self._evaluate_theme(
                aligner=aligner,
                runtime=runtime,
                thematic_ids=thematic_ids,
                theme_id=theme_id,
                original_geometry=feat.geometry,
            )

        return AlignerResult(runtime.process_results_evaluated)

    def _resolve_thematic_ids(
        self,
        *,
        aligner: "Aligner",
        thematic_ids: Optional[List[InputId]],
    ) -> tuple[list[InputId], bool]:
        calculate_zeros = True
        if thematic_ids is None:
            thematic_ids = list(aligner.thematic_data.features.keys())
            calculate_zeros = False
        else:
            thematic_ids = list(thematic_ids)
        return thematic_ids, calculate_zeros

    def _resolve_relevant_distances(
        self, relevant_distances: Optional[Iterable[float]]
    ) -> list[float]:
        if relevant_distances is None:
            return [
                round(k, RELEVANT_DISTANCE_DECIMALS)
                for k in np.arange(0, 310, 10, dtype=int) / 100
            ]
        return list(relevant_distances)

    def _validate_evaluation_inputs(
        self,
        *,
        aligner: "Aligner",
        thematic_ids: list[InputId],
        relevant_distances: list[float],
    ) -> None:
        if any(
            id_to_evaluate not in aligner.thematic_data.features.keys()
            for id_to_evaluate in thematic_ids
        ):
            raise ValueError("not all ids are found in the thematic data")
        if 0 not in relevant_distances:
            raise ValueError(
                "Evaluation cannot be executed when 0 is not available in the array of relevant distances"
            )

    def _prepare_process_results(
        self,
        *,
        aligner: "Aligner",
        thematic_ids: list[InputId],
        relevant_distances: list[float],
        calculate_zeros: bool,
    ) -> tuple[dict, dict]:
        aligner_result = aligner.predict(
            thematic_ids=thematic_ids,
            relevant_distances=relevant_distances,
            diff_metric=aligner.diff_metric,
        )
        process_results = aligner_result.get_results(aligner=aligner)
        process_results_predictions = aligner_result.get_results(
            aligner=aligner, result_type=AlignerResultType.PREDICTIONS
        )

        if calculate_zeros:
            aligner_result_zero = aligner.process(relevant_distances=[0])
            process_results_zero = aligner_result_zero.get_results(aligner=aligner)
            process_results = deep_merge(process_results_zero, process_results)

        return process_results, process_results_predictions

    def _evaluate_theme(
        self,
        *,
        aligner: "Aligner",
        runtime: _EvaluationRuntime,
        thematic_ids: list[InputId],
        theme_id: InputId,
        original_geometry: BaseGeometry,
    ) -> None:
        if theme_id not in thematic_ids:
            runtime.process_results_evaluated = self.update_evaluation_with_original(
                aligner=aligner,
                metadata_field=runtime.metadata_field,
                original_geometry=original_geometry,
                process_results_evaluated=runtime.process_results_evaluated,
                theme_id=theme_id,
                evaluation=Evaluation.NOT_EVALUATED,
            )
            return
        self._evaluate_selected_theme(
            aligner=aligner,
            runtime=runtime,
            theme_id=theme_id,
            original_geometry=original_geometry,
        )

    def _evaluate_selected_theme(
        self,
        *,
        aligner: "Aligner",
        runtime: _EvaluationRuntime,
        theme_id: InputId,
        original_geometry: BaseGeometry,
    ) -> None:
        dict_predictions_results = runtime.process_results_predictions.get(theme_id, {})
        if dict_predictions_results == {}:
            runtime.process_results_evaluated = self.update_evaluation_with_original(
                aligner=aligner,
                metadata_field=runtime.metadata_field,
                original_geometry=original_geometry,
                process_results_evaluated=runtime.process_results_evaluated,
                theme_id=theme_id,
                evaluation=Evaluation.TO_CHECK_NO_PREDICTION,
            )
            return

        if not self._prediction_score_available(
            process_results_temp_predictions=runtime.process_results_temp_predictions,
            process_results_evaluated=runtime.process_results_evaluated,
            theme_id=theme_id,
        ):
            runtime.process_results_evaluated = self.update_evaluation_with_original(
                aligner=aligner,
                metadata_field=runtime.metadata_field,
                original_geometry=original_geometry,
                process_results_evaluated=runtime.process_results_evaluated,
                theme_id=theme_id,
                evaluation=Evaluation.NOT_EVALUATED,
            )
            return

        scores = []
        distances = []
        predictions = []
        observation_match = False
        latest_props = None
        base_brdr_observation = self.get_brdr_observation_from_properties(
            aligner=aligner,
            id_theme=theme_id,
            base_metadata_field=runtime.metadata_field,
        )

        for dist in sorted(dict_predictions_results.keys()):
            outcome = self._score_prediction_candidate(
                aligner=aligner,
                runtime=runtime,
                theme_id=theme_id,
                dist=dist,
                prediction_result=dict_predictions_results[dist],
                base_brdr_observation=base_brdr_observation,
                scores=scores,
                distances=distances,
                predictions=predictions,
            )
            status = outcome["status"]
            latest_props = outcome["props"]
            if status == "observation_match":
                observation_match = True

        best_ix = sorted(range(len(scores)), reverse=True, key=lambda i: scores[i])
        len_best_ix = len(best_ix)
        if not self._handle_no_observation_match(
            runtime=runtime,
            theme_id=theme_id,
            best_ix=best_ix,
            len_best_ix=len_best_ix,
            observation_match=observation_match,
            predictions=predictions,
            latest_props=latest_props,
        ):
            return

        if runtime.max_predictions > 0 and len_best_ix > runtime.max_predictions:
            best_ix = best_ix[: runtime.max_predictions]
        if len(best_ix) > 0:
            for ix in best_ix:
                distance = distances[ix]
                prediction = predictions[ix]
                runtime.process_results_evaluated[theme_id][distance] = prediction
        else:
            runtime.process_results_evaluated = self.update_evaluation_with_original(
                aligner=aligner,
                metadata_field=runtime.metadata_field,
                original_geometry=original_geometry,
                process_results_evaluated=runtime.process_results_evaluated,
                theme_id=theme_id,
                evaluation=Evaluation.TO_CHECK_NO_PREDICTION,
            )

    def _prediction_score_available(
        self,
        *,
        process_results_temp_predictions: dict,
        process_results_evaluated: dict,
        theme_id: InputId,
    ) -> bool:
        for dist in process_results_temp_predictions[theme_id]:
            if not PREDICTION_SCORE in process_results_evaluated[theme_id][dist][
                "properties"
            ]:
                return False
        return True

    def _score_prediction_candidate(
        self,
        *,
        aligner: "Aligner",
        runtime: _EvaluationRuntime,
        theme_id: InputId,
        dist: float,
        prediction_result: dict,
        base_brdr_observation: dict | None,
        scores: list[float],
        distances: list[float],
        predictions: list[dict],
    ) -> dict:
        props = deepcopy(runtime.process_results_evaluated[theme_id][dist]["properties"])
        props_evaluation = self.get_observation_comparison_properties(
            aligner=aligner,
            process_result=prediction_result,
            base_brdr_observation=base_brdr_observation,
        )
        props.update(props_evaluation)

        full = props[FULL_ACTUAL_FIELD_NAME]
        if (
            runtime.full_reference_strategy == FullReferenceStrategy.ONLY_FULL_REFERENCE
            and not full
        ):
            return {"status": "ignored", "props": props}
        if (
            props[EVALUATION_FIELD_NAME]
            in (Evaluation.TO_CHECK_NO_PREDICTION, Evaluation.NOT_EVALUATED)
            and props[PREDICTION_COUNT] == 1
        ):
            props[EVALUATION_FIELD_NAME] = Evaluation.PREDICTION_UNIQUE
        elif (
            props[EVALUATION_FIELD_NAME]
            in (Evaluation.TO_CHECK_NO_PREDICTION, Evaluation.NOT_EVALUATED)
            and props[PREDICTION_COUNT] > 1
        ):
            props[EVALUATION_FIELD_NAME] = Evaluation.TO_CHECK_PREDICTION_MULTI
        elif props[EVALUATION_FIELD_NAME] not in (
            Evaluation.TO_CHECK_NO_PREDICTION,
            Evaluation.NOT_EVALUATED,
        ):
            props[PREDICTION_SCORE] = 100
            scores.clear()
            distances.clear()
            predictions.clear()
            scores.append(props[PREDICTION_SCORE])
            distances.append(dist)
            runtime.process_results_evaluated[theme_id][dist]["properties"] = props
            predictions.append(runtime.process_results_evaluated[theme_id][dist])
            return {"status": "observation_match", "props": props}
        if full:
            if runtime.full_reference_strategy != FullReferenceStrategy.NO_FULL_REFERENCE:
                props[EVALUATION_FIELD_NAME] = Evaluation.TO_CHECK_PREDICTION_FULL
                prediction_score = props[PREDICTION_SCORE] + 50
                if prediction_score > 100:
                    prediction_score = 100
                props[PREDICTION_SCORE] = prediction_score
            else:
                props[EVALUATION_FIELD_NAME] = Evaluation.TO_CHECK_PREDICTION_MULTI_FULL

        scores.append(props[PREDICTION_SCORE])
        distances.append(dist)
        runtime.process_results_temp_predictions[theme_id][dist]["properties"] = props
        predictions.append(runtime.process_results_temp_predictions[theme_id][dist])
        return {"status": "scored", "props": props}

    def _handle_no_observation_match(
        self,
        *,
        runtime: _EvaluationRuntime,
        theme_id: InputId,
        best_ix: list[int],
        len_best_ix: int,
        observation_match: bool,
        predictions: list[dict],
        latest_props: dict | None,
    ) -> bool:
        if observation_match:
            return True

        if len_best_ix == 1:
            props = predictions[0]["properties"]
            if FULL_ACTUAL_FIELD_NAME in props and props[FULL_ACTUAL_FIELD_NAME]:
                predictions[0]["properties"][
                    EVALUATION_FIELD_NAME
                ] = Evaluation.PREDICTION_UNIQUE_AND_FULL_REFERENCE
            else:
                predictions[0]["properties"][EVALUATION_FIELD_NAME] = (
                    Evaluation.PREDICTION_UNIQUE
                )

        if (
            len_best_ix > 1
            and runtime.max_predictions == 1
            and not runtime.multi_to_best_prediction
        ):
            relevant_distance = round(0, RELEVANT_DISTANCE_DECIMALS)
            props = latest_props if latest_props is not None else {}
            props[EVALUATION_FIELD_NAME] = Evaluation.TO_CHECK_ORIGINAL
            props[PREDICTION_SCORE] = -1
            if REMARK_FIELD_NAME in props:
                remarks = props[REMARK_FIELD_NAME]
            else:
                remarks = []
            remarks.append(ProcessRemark.MULTIPLE_PREDICTIONS_ORIGINAL_RETURNED)
            props[REMARK_FIELD_NAME] = remarks
            runtime.process_results_evaluated[theme_id][relevant_distance][
                "properties"
            ].update(props)
            return False

        return True

    def update_evaluation_with_original(
        self,
        *,
        aligner: "Aligner",
        metadata_field: str,
        original_geometry: BaseGeometry,
        process_results_evaluated,
        theme_id: str | int,
        evaluation: Evaluation,
    ) -> Any:
        relevant_distance = round(0, RELEVANT_DISTANCE_DECIMALS)
        try:
            process_result = process_results_evaluated[theme_id][relevant_distance]
            props = deepcopy(process_result["properties"])
        except Exception:
            process_result = {"result": original_geometry}
            props = {}
        base_brdr_observation = self.get_brdr_observation_from_properties(
            aligner=aligner,
            id_theme=theme_id,
            base_metadata_field=metadata_field,
        )
        props_evaluation = self.get_observation_comparison_properties(
            aligner=aligner,
            process_result=process_result,
            base_brdr_observation=base_brdr_observation,
        )
        props.update(props_evaluation)
        props[EVALUATION_FIELD_NAME] = evaluation
        props[PREDICTION_SCORE] = -1
        if REMARK_FIELD_NAME in props:
            remarks = props[REMARK_FIELD_NAME]
        else:
            remarks = []
        remarks.append(ProcessRemark.NOT_EVALUATED_ORIGINAL_RETURNED)
        props[REMARK_FIELD_NAME] = remarks
        process_results_evaluated[theme_id][relevant_distance]["properties"].update(
            props
        )
        return process_results_evaluated

    def get_observation_properties(
        self,
        *,
        aligner: "Aligner",
        process_result,
    ) -> dict:
        geom_process_result = process_result["result"]
        properties = {
            FULL_ACTUAL_FIELD_NAME: None,
        }
        actual_brdr_observation = process_result.get(
            "observations"
        ) or aligner.compare_to_reference(geom_process_result)
        process_result["observations"] = actual_brdr_observation
        if (
            actual_brdr_observation is None
            or geom_process_result is None
            or geom_process_result.is_empty
        ):
            return properties
        properties[FULL_ACTUAL_FIELD_NAME] = actual_brdr_observation["full"]
        return properties

    def get_observation_comparison_properties(
        self,
        *,
        aligner: "Aligner",
        process_result,
        base_brdr_observation=None,
    ) -> dict:
        def _resolve_measure_type(base_obs, actual_obs) -> str:
            for obs in (actual_obs, base_obs):
                if isinstance(obs, dict):
                    measure_type = obs.get("measure_type")
                    if measure_type in {"area", "length", "count"}:
                        return measure_type
            return "area"

        def _metric_value(ref_dict: dict, measure_type: str) -> float:
            if not isinstance(ref_dict, dict):
                return 0.0
            if measure_type in {"area", "length", "count"} and measure_type in ref_dict:
                return float(ref_dict[measure_type])
            per_feature_type = ref_dict.get("measure_type")
            if per_feature_type in {"area", "length", "count"} and per_feature_type in ref_dict:
                return float(ref_dict[per_feature_type])
            for key in ("area", "length", "count"):
                if key in ref_dict:
                    return float(ref_dict[key])
            value = ref_dict.get(measure_type)
            if value is None:
                value = ref_dict.get("area", 0.0)
            return float(value)

        def _pair_metric_values(base_ref: dict, actual_ref: dict, default_type: str):
            candidate_order = [default_type, "area", "length", "count"]
            seen = set()
            ordered = []
            for c in candidate_order:
                if c in {"area", "length", "count"} and c not in seen:
                    ordered.append(c)
                    seen.add(c)
            for key in ordered:
                if key in base_ref and key in actual_ref:
                    return float(base_ref[key]), float(actual_ref[key])
            return _metric_value(base_ref, default_type), _metric_value(
                actual_ref, default_type
            )

        def _od_value(obs: dict, measure_type: str) -> float | None:
            od = obs.get("reference_od")
            if od is None:
                return None
            if not isinstance(od, dict):
                return None
            if measure_type in od:
                return float(od[measure_type])
            for key in ("area", "length", "count"):
                if key in od:
                    return float(od[key])
            value = od.get(measure_type)
            if value is None:
                value = od.get("area")
            return None if value is None else float(value)

        geom_process_result = process_result["result"]
        threshold_od_percentage = 1
        properties = {
            EVALUATION_FIELD_NAME: Evaluation.TO_CHECK_NO_PREDICTION,
            FULL_BASE_FIELD_NAME: None,
            FULL_ACTUAL_FIELD_NAME: None,
            OD_ALIKE_FIELD_NAME: None,
            EQUAL_REFERENCE_FEATURES_FIELD_NAME: None,
            DIFF_PERCENTAGE_FIELD_NAME: None,
            DIFF_AREA_FIELD_NAME: None,
        }
        props = self.get_observation_properties(
            aligner=aligner, process_result=process_result
        )
        properties.update(props)

        actual_brdr_observation = process_result.get(
            "observations"
        ) or aligner.compare_to_reference(geom_process_result)
        if (
            actual_brdr_observation is None
            or geom_process_result is None
            or geom_process_result.is_empty
        ):
            return properties

        if not is_brdr_observation(base_brdr_observation):
            properties[EVALUATION_FIELD_NAME] = Evaluation.TO_CHECK_NO_PREDICTION
            return properties
        properties[FULL_BASE_FIELD_NAME] = base_brdr_observation["full"]
        measure_type = _resolve_measure_type(
            base_obs=base_brdr_observation,
            actual_obs=actual_brdr_observation,
        )
        od_alike = False
        base_od = _od_value(base_brdr_observation, measure_type)
        actual_od = _od_value(actual_brdr_observation, measure_type)
        if base_od is None and actual_od is None:
            od_alike = True
        elif base_od is None or actual_od is None:
            od_alike = False
        elif base_od == 0:
            od_alike = actual_od == 0
        elif (abs(base_od - actual_od) * 100 / base_od) < threshold_od_percentage:
            od_alike = True
        properties[OD_ALIKE_FIELD_NAME] = od_alike

        equal_reference_features = False
        if (
            base_brdr_observation["reference_features"].keys()
            == actual_brdr_observation["reference_features"].keys()
        ):
            equal_reference_features = True
            max_diff_area_reference_feature = 0
            max_diff_percentage_reference_feature = 0
            for key in base_brdr_observation["reference_features"].keys():
                if (
                    base_brdr_observation["reference_features"][key]["full"]
                    != actual_brdr_observation["reference_features"][key]["full"]
                ):
                    equal_reference_features = False

                base_measure, actual_measure = _pair_metric_values(
                    base_brdr_observation["reference_features"][key],
                    actual_brdr_observation["reference_features"][key],
                    measure_type,
                )
                diff_area_reference_feature = abs(
                    base_measure - actual_measure
                )
                area = base_measure
                if area > 0:
                    diff_percentage_reference_feature = (
                        abs(base_measure - actual_measure) * 100 / area
                    )
                else:
                    diff_percentage_reference_feature = 0
                if diff_area_reference_feature > max_diff_area_reference_feature:
                    max_diff_area_reference_feature = diff_area_reference_feature
                if (
                    diff_percentage_reference_feature
                    > max_diff_percentage_reference_feature
                ):
                    max_diff_percentage_reference_feature = (
                        diff_percentage_reference_feature
                    )
            properties[EQUAL_REFERENCE_FEATURES_FIELD_NAME] = equal_reference_features
            properties[DIFF_AREA_FIELD_NAME] = max_diff_area_reference_feature
            properties[DIFF_PERCENTAGE_FIELD_NAME] = (
                max_diff_percentage_reference_feature
            )
        # EVALUATION
        if (
            equal_reference_features
            and od_alike
            and base_brdr_observation["full"]
            and actual_brdr_observation["full"]
        ):  # observation is the same, and both geometries are 'full'
            properties[EVALUATION_FIELD_NAME] = (
                Evaluation.EQUALITY_BY_ID_AND_FULL_REFERENCE
            )
        elif (
            equal_reference_features
            and od_alike
            and base_brdr_observation["full"] == actual_brdr_observation["full"]
        ):  # observation is the same,  both geometries are not 'full'
            properties[EVALUATION_FIELD_NAME] = Evaluation.EQUALITY_BY_ID
        elif (
            base_brdr_observation["full"]
            and actual_brdr_observation["full"]
            and od_alike
        ):  # observation not the same but geometries are full
            properties[EVALUATION_FIELD_NAME] = Evaluation.EQUALITY_BY_FULL_REFERENCE
        else:
            properties[EVALUATION_FIELD_NAME] = Evaluation.TO_CHECK_NO_PREDICTION
        return properties

    def get_brdr_observation_from_properties(
        self,
        *,
        aligner: "Aligner",
        id_theme: Any,
        base_metadata_field: str,
    ) -> dict:
        try:
            base_metadata = aligner.thematic_data.features.get(id_theme).properties[
                base_metadata_field
            ]
            if isinstance(base_metadata, str):
                base_metadata = json.loads(base_metadata)

            base_brdr_observation = _reverse_metadata_observations_to_brdr_observation(
                base_metadata
            )
        except Exception:
            base_brdr_observation = None
        return base_brdr_observation
