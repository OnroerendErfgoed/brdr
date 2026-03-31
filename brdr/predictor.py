from abc import ABC
from abc import abstractmethod
from typing import List
from typing import Optional
from typing import TYPE_CHECKING
from typing import Union

import numpy as np

from brdr.constants import PREDICTION_COUNT
from brdr.constants import PREDICTION_SCORE
from brdr.constants import RELEVANT_DISTANCE_DECIMALS
from brdr.constants import STABILITY
from brdr.constants import ZERO_STREAK
from brdr.enums import DiffMetric
from brdr.logger import Logger
from brdr.typings import InputId
from brdr.utils import coverage_ratio
from brdr.utils import determine_stability
from brdr.utils import get_geometry_difference_metrics_from_processresults

if TYPE_CHECKING:
    from brdr.aligner import Aligner
    from brdr.aligner import AlignerResult


class BasePredictor(ABC):
    """
    Abstract base class for prediction strategies.

    A predictor analyzes process results across a series of relevant distances
    and marks interesting candidate results with prediction metadata.
    """

    def __init__(self, feedback=None):
        self.logger = Logger(feedback)

    @abstractmethod
    def predict(
        self,
        *,
        aligner: "Aligner",
        relevant_distances: Optional[Union[List[float], np.ndarray]] = None,
        thematic_ids: Optional[List[InputId]] = None,
        diff_metric: Optional[DiffMetric] = None,
    ) -> "AlignerResult":
        """Execute prediction strategy for the given aligner context."""
        pass


class AlignerPredictor(BasePredictor):
    """
    Default predictor implementation with the original Aligner prediction logic.
    """

    def predict(
        self,
        *,
        aligner: "Aligner",
        relevant_distances: Optional[Union[List[float], np.ndarray]] = None,
        thematic_ids: Optional[List[InputId]] = None,
        diff_metric: Optional[DiffMetric] = None,
    ) -> "AlignerResult":
        from brdr.aligner import AlignerResult

        if thematic_ids is None:
            thematic_ids = aligner.thematic_data.features.keys()
        if any(
            id_to_predict not in aligner.thematic_data.features.keys()
            for id_to_predict in thematic_ids
        ):
            raise ValueError("not all ids are found in the thematic data")
        if relevant_distances is None:
            relevant_distances = [
                round(k, RELEVANT_DISTANCE_DECIMALS)
                for k in np.arange(0, 310, 10, dtype=int) / 100
            ]
        rd_prediction = list(relevant_distances)
        max_relevant_distance = max(rd_prediction)
        cvg_ratio = coverage_ratio(values=relevant_distances, min_val=0, bin_count=10)
        cvg_ratio_threshold = 0.75
        # cvg_ratio: indication of the rd values can be used to make a brdr_prediction_score. When there is enough coverage of predictions to determine a prediction_score we also add 0 and a long-range value(+1).
        # Otherwise we only add a short-range value (+0.1) to check for stability
        if cvg_ratio > cvg_ratio_threshold:
            rd_prediction.append(round(0, RELEVANT_DISTANCE_DECIMALS))
            rd_prediction.append(
                round(max_relevant_distance + 0.1, RELEVANT_DISTANCE_DECIMALS)
            )
            rd_prediction.append(
                round(max_relevant_distance + 1, RELEVANT_DISTANCE_DECIMALS)
            )
        else:
            rd_prediction.append(
                round(max_relevant_distance + 0.1, RELEVANT_DISTANCE_DECIMALS)
            )
        rd_prediction = list(set(rd_prediction))
        rd_prediction = sorted(rd_prediction)

        # Get aligner_result for all relevant_distances
        aligner_result = aligner.process(
            thematic_ids=thematic_ids,
            relevant_distances=rd_prediction,
        )
        process_results = aligner_result.results

        # Search for predictions
        if diff_metric is None:
            diff_metric = aligner.diff_metric
        for theme_id, process_result in process_results.items():
            diffs = get_geometry_difference_metrics_from_processresults(
                process_result,
                aligner.thematic_data.features.get(theme_id).geometry,
                aligner.reference_data.union,
                diff_metric=diff_metric,
            )
            if len(diffs) != len(rd_prediction):
                self.logger.feedback_warning(
                    f"Number of computed diffs for thematic element {theme_id} does "
                    f"not match the number of relevant distances."
                )
                continue
            diff_values = list(diffs.values())
            dict_stability = determine_stability(rd_prediction, diff_values)
            prediction_count = 0
            for rd in rd_prediction:
                if rd not in relevant_distances:
                    del process_results[theme_id][rd]
                    continue

                process_results[theme_id][rd]["properties"][STABILITY] = dict_stability[
                    rd
                ][STABILITY]
                if dict_stability[rd][ZERO_STREAK] is not None:
                    if cvg_ratio > cvg_ratio_threshold:
                        prediction_count += 1
                        process_results[theme_id][rd]["properties"][
                            PREDICTION_SCORE
                        ] = dict_stability[rd][ZERO_STREAK][3]
            for rd, process_result in process_results[theme_id].items():
                if PREDICTION_SCORE in process_result["properties"]:
                    process_result["properties"][PREDICTION_COUNT] = prediction_count
        return AlignerResult(process_results)
