import numpy as np

from brdr.aligner import Aligner, AlignerResult
from brdr.constants import EVALUATION_FIELD_NAME
from brdr.constants import RELEVANT_DISTANCE_FIELD_NAME
from brdr.evaluator import AlignerEvaluator
from brdr.loader import DictLoader
from shapely import from_wkt


class TaggedEvaluator(AlignerEvaluator):
    """
    Add a custom tag to each evaluated process result.
    """

    def evaluate(self, **kwargs) -> AlignerResult:
        result = super().evaluate(**kwargs)
        for _, by_distance in result.results.items():
            for _, process_result in by_distance.items():
                if process_result is None:
                    continue
                process_result["properties"]["custom_evaluator_tag"] = "checked"
        return result


if __name__ == "__main__":
    """Example: inject a custom evaluator into Aligner."""

    thematic = {
        "theme_1": from_wkt("POLYGON ((0 0, 0 9, 5 10, 10 0, 0 0))"),
    }
    reference = {
        "ref_1": from_wkt("POLYGON ((0 1, 0 10, 8 10, 10 1, 0 1))"),
    }

    aligner = Aligner(evaluator=TaggedEvaluator())
    aligner.load_thematic_data(DictLoader(thematic))
    aligner.load_reference_data(DictLoader(reference))

    distances = np.arange(0, 210, 10, dtype=int) / 100
    evaluated = aligner.evaluate(
        relevant_distances=distances,
        metadata_field=None,
    ).get_results(aligner=aligner)

    for _, by_distance in evaluated.items():
        for _, process_result in by_distance.items():
            if process_result is None:
                continue
            print(
                f"distance={process_result['properties'][RELEVANT_DISTANCE_FIELD_NAME]}, "
                f"evaluation={process_result['properties'].get(EVALUATION_FIELD_NAME)}, "
                f"tag={process_result['properties'].get('custom_evaluator_tag')}"
            )
