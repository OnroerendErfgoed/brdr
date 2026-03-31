import numpy as np

from brdr.aligner import Aligner, AlignerResult
from brdr.constants import PREDICTION_SCORE
from brdr.loader import DictLoader
from brdr.predictor import AlignerPredictor
from shapely import from_wkt


class TopPredictionOnlyPredictor(AlignerPredictor):
    """
    Keep only the highest-scoring prediction per thematic feature.
    """

    def predict(self, **kwargs) -> AlignerResult:
        result = super().predict(**kwargs)
        for theme_id, by_distance in result.results.items():
            with_score = [
                (rd, pr)
                for rd, pr in by_distance.items()
                if pr and PREDICTION_SCORE in pr["properties"]
            ]
            if not with_score:
                continue
            best_rd, _ = max(with_score, key=lambda item: item[1]["properties"][PREDICTION_SCORE])
            for rd in list(by_distance.keys()):
                if rd != best_rd:
                    del by_distance[rd]
        return result


if __name__ == "__main__":
    """Example: inject a custom predictor into Aligner."""

    thematic = {
        "theme_1": from_wkt("POLYGON ((0 0, 0 9, 5 10, 10 0, 0 0))"),
    }
    reference = {
        "ref_1": from_wkt("POLYGON ((0 1, 0 10, 8 10, 10 1, 0 1))"),
    }

    aligner = Aligner(predictor=TopPredictionOnlyPredictor())
    aligner.load_thematic_data(DictLoader(thematic))
    aligner.load_reference_data(DictLoader(reference))

    distances = np.arange(0, 210, 10, dtype=int) / 100
    predictions = aligner.predict(relevant_distances=distances).results["theme_1"]

    print("retained distances:", sorted(predictions.keys()))
    for rd, process_result in predictions.items():
        score = process_result["properties"].get(PREDICTION_SCORE)
        print(f"distance={rd}, prediction_score={score}")
