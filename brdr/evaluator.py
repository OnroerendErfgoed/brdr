import json
from abc import ABC
from abc import abstractmethod
from copy import deepcopy
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
from brdr.typings import InputId
from brdr.utils import deep_merge
from brdr.utils import is_brdr_observation

if TYPE_CHECKING:
    from brdr.aligner import Aligner
    from brdr.aligner import AlignerResult


def _reverse_metadata_observations_to_brdr_observation(metadata: dict) -> dict:
    observations = metadata["observations"]
    result_uri = observations[0]["used"] if observations else None

    brdr_observation = {
        "area": None,
        "full": None,
        "reference_features": {},
        "reference_od": None,
        "result": result_uri,
    }

    for obs in observations:
        prop = obs.get("observed_property", "")
        result = obs.get("result", {}).get("value")
        target = obs.get("has_feature_of_interest")

        if prop.endswith("area_overlap"):
            ref = brdr_observation["reference_features"].setdefault(target, {})
            ref["area"] = result
        elif prop.endswith("area_overlap_percentage"):
            ref = brdr_observation["reference_features"].setdefault(target, {})
            ref["percentage"] = result
        elif prop.endswith("area_overlap_full"):
            brdr_observation["full"] = result
        elif prop.endswith("area_open_domain"):
            brdr_observation["reference_od"] = {"area": result}
        elif prop.endswith("area"):
            brdr_observation["area"] = result

    return brdr_observation


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

        calculate_zeros = True  # boolean to check if we need to add al zero-rd results to the evaluations
        if thematic_ids is None:
            thematic_ids = aligner.thematic_data.features.keys()
            calculate_zeros = False  # when all thematic features will be calculated, there is no need to calculate the zeros for all seperately
        if any(
            id_to_evaluate not in aligner.thematic_data.features.keys()
            for id_to_evaluate in thematic_ids
        ):
            raise ValueError("not all ids are found in the thematic data")
        if relevant_distances is None:
            relevant_distances = [
                round(k, RELEVANT_DISTANCE_DECIMALS)
                for k in np.arange(0, 310, 10, dtype=int) / 100
            ]

        if 0 not in relevant_distances:
            raise ValueError(
                "Evaluation cannot be executed when 0 is not available in the array of relevant distances"
            )

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
            # Calculate the ZERO-situation for all, only when the predict is not done for all thematic ids
            aligner_result_zero = aligner.process(
                relevant_distances=[0],
            )
            process_results_zero = aligner_result_zero.get_results(aligner=aligner)
            process_results = deep_merge(
                process_results_zero,
                process_results,
            )
        process_results_temp_predictions = deepcopy(process_results_predictions)
        process_results_evaluated = deepcopy(process_results)

        for theme_id, feat in aligner.thematic_data.features.items():
            original_geometry = feat.geometry
            if theme_id not in thematic_ids:
                # PART 1)NOT EVALUATED
                process_results_evaluated = self.update_evaluation_with_original(
                    aligner=aligner,
                    metadata_field=metadata_field,
                    original_geometry=original_geometry,
                    process_results_evaluated=process_results_evaluated,
                    theme_id=theme_id,
                    evaluation=Evaluation.NOT_EVALUATED,
                )

            # Features are split up in 2 groups: TO_EVALUATE and NOT_TO_EVALUATE (original returned)
            # The evaluated features will be split up:
            #   *No prediction available
            #   *Predictions available

            # PART 2: TO_EVALUATE
            elif theme_id in thematic_ids:

                if (
                    theme_id not in process_results_predictions.keys()
                    or process_results_predictions[theme_id] == {}
                ):
                    # No predictions available
                    process_results_evaluated = self.update_evaluation_with_original(
                        aligner=aligner,
                        metadata_field=metadata_field,
                        original_geometry=original_geometry,
                        process_results_evaluated=process_results_evaluated,
                        theme_id=theme_id,
                        evaluation=Evaluation.TO_CHECK_NO_PREDICTION,
                    )
                    continue
                # Check if prediction_scores are available to do the evaluation
                prediction_score_available = True
                for dist in process_results_temp_predictions[theme_id]:
                    if (
                        not PREDICTION_SCORE
                        in process_results_evaluated[theme_id][dist]["properties"]
                    ):
                        prediction_score_available = False
                # thematic objects that do not have a prediction_score from predict() are not evaluated and returned as they are
                if not prediction_score_available:
                    process_results_evaluated = self.update_evaluation_with_original(
                        aligner=aligner,
                        metadata_field=metadata_field,
                        original_geometry=original_geometry,
                        process_results_evaluated=process_results_evaluated,
                        theme_id=theme_id,
                        evaluation=Evaluation.NOT_EVALUATED,
                    )
                    continue

                # When there are predictions available
                dict_predictions_results = process_results_predictions[theme_id]
                scores = []
                distances = []
                predictions = []
                observation_match = False
                base_brdr_observation = self.get_brdr_observation_from_properties(
                    aligner=aligner,
                    id_theme=theme_id,
                    base_metadata_field=metadata_field,
                )
                for dist in sorted(dict_predictions_results.keys()):
                    props = deepcopy(
                        process_results_evaluated[theme_id][dist]["properties"]
                    )
                    props_evaluation = self.get_observation_comparison_properties(
                        aligner=aligner,
                        process_result=dict_predictions_results[dist],
                        base_brdr_observation=base_brdr_observation,
                    )
                    props.update(props_evaluation)

                    full = props[FULL_ACTUAL_FIELD_NAME]
                    if (
                        full_reference_strategy
                        == FullReferenceStrategy.ONLY_FULL_REFERENCE
                        and not full
                    ):
                        # this prediction is ignored
                        continue
                    if (
                        props[EVALUATION_FIELD_NAME]
                        in (Evaluation.TO_CHECK_NO_PREDICTION, Evaluation.NOT_EVALUATED)
                        and props[PREDICTION_COUNT] == 1
                    ):
                        props[EVALUATION_FIELD_NAME] = Evaluation.PREDICTION_UNIQUE
                        # TODO can we add continue here? No, because this can get overwritten by prediction_unique_full
                    elif (
                        props[EVALUATION_FIELD_NAME]
                        in (Evaluation.TO_CHECK_NO_PREDICTION, Evaluation.NOT_EVALUATED)
                        and props[PREDICTION_COUNT] > 1
                    ):
                        props[EVALUATION_FIELD_NAME] = (
                            Evaluation.TO_CHECK_PREDICTION_MULTI
                        )
                    elif props[EVALUATION_FIELD_NAME] not in (
                        Evaluation.TO_CHECK_NO_PREDICTION,
                        Evaluation.NOT_EVALUATED,
                    ):
                        # this prediction has a equality based on observation so the rest is not checked anymore
                        observation_match = True
                        props[PREDICTION_SCORE] = 100
                        scores = []
                        distances = []
                        predictions = []
                        scores.append(props[PREDICTION_SCORE])
                        distances.append(dist)
                        process_results_evaluated[theme_id][dist]["properties"] = props
                        predictions.append(process_results_evaluated[theme_id][dist])
                        continue
                    if full:
                        if (
                            full_reference_strategy
                            != FullReferenceStrategy.NO_FULL_REFERENCE
                        ):
                            props[EVALUATION_FIELD_NAME] = (
                                Evaluation.TO_CHECK_PREDICTION_FULL
                            )
                            prediction_score = props[PREDICTION_SCORE] + 50
                            if prediction_score > 100:
                                prediction_score = 100
                            props[PREDICTION_SCORE] = prediction_score
                        else:
                            props[EVALUATION_FIELD_NAME] = (
                                Evaluation.TO_CHECK_PREDICTION_MULTI_FULL
                            )

                    scores.append(props[PREDICTION_SCORE])
                    distances.append(dist)
                    process_results_temp_predictions[theme_id][dist][
                        "properties"
                    ] = props
                    predictions.append(process_results_temp_predictions[theme_id][dist])

                # get max amount of best-scoring predictions
                best_ix = sorted(
                    range(len(scores)), reverse=True, key=lambda i: scores[i]
                )
                len_best_ix = len(best_ix)

                if not observation_match:
                    # if there is only one prediction left,  evaluation is set to PREDICTION_UNIQUE_FULL
                    if len_best_ix == 1 and not observation_match:
                        props = predictions[0]["properties"]
                        if (
                            FULL_ACTUAL_FIELD_NAME in props
                            and props[FULL_ACTUAL_FIELD_NAME]
                        ):
                            predictions[0]["properties"][
                                EVALUATION_FIELD_NAME
                            ] = Evaluation.PREDICTION_UNIQUE_AND_FULL_REFERENCE
                        else:
                            predictions[0]["properties"][
                                EVALUATION_FIELD_NAME
                            ] = Evaluation.PREDICTION_UNIQUE

                    # if there are multiple predictions, but we want only one and we ask for the original
                    if (
                        len_best_ix > 1
                        and max_predictions == 1
                        and not multi_to_best_prediction
                        and not observation_match
                    ):
                        relevant_distance = round(0, RELEVANT_DISTANCE_DECIMALS)
                        props[EVALUATION_FIELD_NAME] = Evaluation.TO_CHECK_ORIGINAL
                        props[PREDICTION_SCORE] = -1
                        if REMARK_FIELD_NAME in props:
                            remarks = props[REMARK_FIELD_NAME]
                        else:
                            remarks = []
                        remarks.append(
                            ProcessRemark.MULTIPLE_PREDICTIONS_ORIGINAL_RETURNED
                        )
                        props[REMARK_FIELD_NAME] = remarks
                        process_results_evaluated[theme_id][relevant_distance][
                            "properties"
                        ].update(props)
                        continue

                if max_predictions > 0 and len_best_ix > max_predictions:
                    best_ix = best_ix[:max_predictions]
                if len(best_ix) > 0:
                    for ix in best_ix:
                        distance = distances[ix]
                        prediction = predictions[ix]
                        process_results_evaluated[theme_id][distance] = prediction
                else:
                    # #when no evaluated predictions, the original is returned
                    process_results_evaluated = self.update_evaluation_with_original(
                        aligner=aligner,
                        metadata_field=metadata_field,
                        original_geometry=original_geometry,
                        process_results_evaluated=process_results_evaluated,
                        theme_id=theme_id,
                        evaluation=Evaluation.TO_CHECK_NO_PREDICTION,
                    )

        return AlignerResult(process_results_evaluated)

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
        od_alike = False
        if (
            base_brdr_observation["reference_od"] is None
            and actual_brdr_observation["reference_od"] is None
        ):
            od_alike = True
        elif (
            base_brdr_observation["reference_od"] is None
            or actual_brdr_observation["reference_od"] is None
        ):
            od_alike = False
        elif (
            abs(
                base_brdr_observation["reference_od"]["area"]
                - actual_brdr_observation["reference_od"]["area"]
            )
            * 100
            / base_brdr_observation["reference_od"]["area"]
        ) < threshold_od_percentage:
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

                diff_area_reference_feature = abs(
                    base_brdr_observation["reference_features"][key]["area"]
                    - actual_brdr_observation["reference_features"][key]["area"]
                )
                area = base_brdr_observation["reference_features"][key]["area"]
                if area > 0:
                    diff_percentage_reference_feature = (
                        abs(
                            base_brdr_observation["reference_features"][key]["area"]
                            - actual_brdr_observation["reference_features"][key]["area"]
                        )
                        * 100
                        / base_brdr_observation["reference_features"][key]["area"]
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
