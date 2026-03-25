from enum import Enum


class OpenDomainStrategy(Enum):
    """
    Strategies for thematic areas that are not covered by the reference layer.

    Attributes
    ----------
    EXCLUDE : str
        Exclude open-domain parts from the resulting geometry.
    AS_IS : str
        Keep open-domain parts unchanged.
    SNAP_INNER_SIDE : str
        Snap open-domain parts to inner reference boundaries.
    SNAP_ALL_SIDE : str
        Snap open-domain parts to both inner and outer sides where applicable.
    SNAP_PREFER_VERTICES : str
        Snap with preference for reference vertices.
    SNAP_NO_PREFERENCE : str
        Snap without vertex/edge preference.
    SNAP_ONLY_VERTICES : str
        Restrict snapping targets to reference vertices only.
    """

    EXCLUDE = "EXCLUDE"
    AS_IS = "AS_IS"
    SNAP_INNER_SIDE = "SNAP_INNER_SIDE"
    SNAP_ALL_SIDE = "SNAP_ALL_SIDE"
    SNAP_PREFER_VERTICES = "SNAP_PREFER_VERTICES"
    SNAP_NO_PREFERENCE = "SNAP_NO_PREFERENCE"
    SNAP_ONLY_VERTICES = "SNAP_ONLY_VERTICES"


class AlignerResultType(str, Enum):
    """
    Output format of alignment runs.

    - `PREDICTIONS`: only predicted candidate results
    - `EVALUATED_PREDICTIONS`: evaluated/scored predictions
    - `PROCESSRESULTS`: full result set for all processed distances

    Attributes
    ----------
    PREDICTIONS : str
        Only predicted candidate results.
    EVALUATED_PREDICTIONS : str
        Predictions with evaluation labels and scores.
    PROCESSRESULTS : str
        Complete processing output for all requested distances.
    """

    PREDICTIONS = "predictions"
    EVALUATED_PREDICTIONS = "evaluated_predictions"
    PROCESSRESULTS = "processresults"


class AlignerInputType(str, Enum):
    """
    Role of an input dataset in the aligner pipeline.

    Attributes
    ----------
    THEMATIC : str
        The input dataset that is being aligned.
    REFERENCE : str
        The target dataset used as alignment reference.
    """

    THEMATIC = "thematic"
    REFERENCE = "reference"


class DiffMetric(str, Enum):
    """
    Metrics used to quantify geometric differences.

    Attributes
    ----------
    AREA_CHANGE : str
        Absolute area difference.
    AREA_PERCENTAGE_CHANGE : str
        Relative area difference in percent.
    SYMMETRICAL_AREA_CHANGE : str
        Symmetric-difference area.
    SYMMETRICAL_AREA_PERCENTAGE_CHANGE : str
        Symmetric-difference area in percent.
    LENGTH_CHANGE : str
        Absolute length difference.
    LENGTH_PERCENTAGE_CHANGE : str
        Relative length difference in percent.
    LENGTH_ADDED_AND_REMOVED : str
        Combined added and removed length.
    LENGTH_REMOVED : str
        Removed length only.
    TOTAL_DISTANCE : str
        Total displacement distance.
    REFERENCE_USAGE : str
        Amount of reference geometry used.
    """

    AREA_CHANGE = "AREA_CHANGE"
    AREA_PERCENTAGE_CHANGE = "AREA_PERCENTAGE_CHANGE"
    SYMMETRICAL_AREA_CHANGE = "SYMMETRICAL_AREA_CHANGE"
    SYMMETRICAL_AREA_PERCENTAGE_CHANGE = "SYMMETRICAL_AREA_PERCENTAGE_CHANGE"
    LENGTH_CHANGE = "LENGTH_CHANGE"
    LENGTH_PERCENTAGE_CHANGE = "LENGTH_PERCENTAGE_CHANGE"
    LENGTH_ADDED_AND_REMOVED = "LENGTH_ADDED_AND_REMOVED"
    LENGTH_REMOVED = "LENGTH_REMOVED"
    TOTAL_DISTANCE = "TOTAL_DISTANCE"
    REFERENCE_USAGE = "REFERENCE_USAGE"


class ProcessRemark(str, Enum):
    """
    Processing status remarks for auditing and diagnostics.

    Attributes
    ----------
    RESULT_UNCHANGED : str
        Result equals original geometry within tolerance.
    INPUT_CIRCLE : str
        Circle-like input detected; original geometry returned.
    RESULT_EMPTY_RETURNED : str
        Result became empty and empty geometry was returned.
    CHANGED_GEOMETRYTYPE_EMPTY_RETURNED : str
        Geometry type changed unexpectedly; empty geometry returned.
    CHANGED_AMOUNT_GEOMETRIES : str
        Number of sub-geometries changed.
    NO_PREDICTION_ORIGINAL_RETURNED : str
        No prediction accepted; original geometry returned.
    MULTIPLE_PREDICTIONS_ORIGINAL_RETURNED : str
        Multiple predictions prevented single choice; original returned.
    NOT_EVALUATED_ORIGINAL_RETURNED : str
        Evaluation not possible; original returned.
    """

    RESULT_UNCHANGED = "resulting geometry equal to original geometry"
    INPUT_CIRCLE = "circle detected: original_geometry_returned"
    RESULT_EMPTY_RETURNED = "resulting geometry empty: empty geometry returned"
    CHANGED_GEOMETRYTYPE_EMPTY_RETURNED = (
        "resulting geometry has different geometrytype: empty geometry returned"
    )
    CHANGED_AMOUNT_GEOMETRIES = (
        "resulting (multi) geometry has different amount of geometries"
    )
    NO_PREDICTION_ORIGINAL_RETURNED = "NO_PREDICTION_ORIGINAL_RETURNED"
    MULTIPLE_PREDICTIONS_ORIGINAL_RETURNED = "MULTIPLE_PREDICTIONS_ORIGINAL_RETURNED"
    NOT_EVALUATED_ORIGINAL_RETURNED = "NOT_EVALUATED_ORIGINAL_RETURNED"


class Evaluation(str, Enum):
    """
    Classification labels for prediction quality and confidence.

    Attributes
    ----------
    EQUALITY_BY_ID_AND_FULL_REFERENCE : str
        Match by reference identifiers and full-reference condition.
    EQUALITY_BY_ID : str
        Match by reference identifiers.
    EQUALITY_BY_FULL_REFERENCE : str
        Match by full-reference condition.
    PREDICTION_UNIQUE : str
        One unique prediction selected.
    PREDICTION_UNIQUE_AND_FULL_REFERENCE : str
        One unique prediction selected with full-reference condition.
    TO_CHECK_PREDICTION_FULL : str
        Manual check advised; prediction selected with full-reference context.
    TO_CHECK_PREDICTION_MULTI : str
        Manual check advised; multiple prediction candidates.
    TO_CHECK_PREDICTION_MULTI_FULL : str
        Manual check advised; multiple candidates with full-reference context.
    TO_CHECK_ORIGINAL : str
        Manual check advised; keep or inspect original geometry.
    TO_CHECK_NO_PREDICTION : str
        Manual check advised; no suitable prediction found.
    NOT_EVALUATED : str
        No evaluation performed.
    NO_CHANGE : str
        No meaningful change detected.
    """

    EQUALITY_BY_ID_AND_FULL_REFERENCE = "equality_by_id_and_full_reference"
    EQUALITY_BY_ID = "equality_by_id"
    EQUALITY_BY_FULL_REFERENCE = "equality_by_full_reference"
    PREDICTION_UNIQUE = "prediction_unique"
    PREDICTION_UNIQUE_AND_FULL_REFERENCE = "prediction_unique_full"
    TO_CHECK_PREDICTION_FULL = "to_check_prediction_full"
    TO_CHECK_PREDICTION_MULTI = "to_check_prediction_multi"
    TO_CHECK_PREDICTION_MULTI_FULL = "to_check_prediction_multi_full"
    TO_CHECK_ORIGINAL = "to_check_original"
    TO_CHECK_NO_PREDICTION = "to_check_no_prediction"
    NOT_EVALUATED = "not_evaluated"
    NO_CHANGE = "no_change"


class ProcessorID(str, Enum):
    """
    Stable identifiers for available processing algorithms.

    Attributes
    ----------
    DIEUSSAERT : str
        Identifier for DieussaertGeometryProcessor.
    SNAP : str
        Identifier for SnapGeometryProcessor.
    NETWORK : str
        Identifier for NetworkGeometryProcessor.
    ALIGNER : str
        Identifier for AlignerGeometryProcessor.
    TOPOLOGY : str
        Identifier for TopologyProcessor.
    """

    DIEUSSAERT = "2024:dieussaert2024a"
    SNAP = "2024:snap2024a"
    NETWORK = "2024:network2024a"
    ALIGNER = "2024:aligner2024a"
    TOPOLOGY = "2024:topology2024a"


class SnapStrategy(str, Enum):
    """
    Priority behavior used during snapping.

    Attributes
    ----------
    ONLY_VERTICES : str
        Snap only to reference vertices.
    PREFER_VERTICES : str
        Prefer vertices over edges.
    NO_PREFERENCE : str
        Use nearest valid target without preference.
    """

    ONLY_VERTICES = "only_vertices"
    PREFER_VERTICES = "prefer_vertices"
    NO_PREFERENCE = "no_preference"


class FullReferenceStrategy(str, Enum):
    """
    How strongly full-reference matches are preferred.

    Attributes
    ----------
    ONLY_FULL_REFERENCE : str
        Keep only full-reference matches.
    PREFER_FULL_REFERENCE : str
        Prefer full-reference matches when available.
    NO_FULL_REFERENCE : str
        Do not prioritize full-reference matches.
    """

    ONLY_FULL_REFERENCE = "only_full_reference"
    PREFER_FULL_REFERENCE = "prefer_full_reference"
    NO_FULL_REFERENCE = "no_full_reference"


class PredictionStrategy(str, Enum):
    """
    How a primary prediction is selected from candidates.

    Attributes
    ----------
    ALL : str
        Keep all prediction candidates.
    BEST : str
        Keep only the best candidate.
    ORIGINAL : str
        Keep the original geometry as final choice.
    """

    ALL = "all"
    BEST = "best"
    ORIGINAL = "original"
