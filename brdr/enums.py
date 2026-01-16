from enum import Enum


class OpenDomainStrategy(Enum):
    """
    Strategies for processing thematic areas not covered by reference data.

    In GIS terms, this defines how the algorithm handles the 'Public Domain'
    when the reference layer (e.g., land parcels) does not provide a target.

    Attributes
    ----------
    EXCLUDE : str
        Completely remove parts of the thematic geometry not on the reference layer.
    AS_IS : str
        Keep parts not covered by reference data unchanged in the result.
    SNAP_INNER_SIDE : str
        Snap everything within the relevant distance to the interior parcel boundary.
    SNAP_ALL_SIDE : str
        Snap to both inner and outer boundaries where possible.
    SNAP_PREFER_VERTICES : str
        Snap to reference polygons, prioritizing vertices over edges.
    SNAP_NO_PREFERENCE : str
        Snap to the closest reference component (edge or vertex) without bias.
    SNAP_ONLY_VERTICES : str
        Only allow snapping to existing vertices of the reference polygons.

    Notes
    -----

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
    2 Types of resulting dictionary are available in Aligner

    * PREDICTIONS = "predictions" (only the predicted versions for specific relevenat distances are returned)
    * PROCESSRESULTS = "processresults" (All versions of resulting geometries for all relevant distances are returned)
    """

    PREDICTIONS = "predictions"
    EVALUATED_PREDICTIONS = "evaluated_predictions"
    PROCESSRESULTS = "processresults"


class AlignerInputType(str, Enum):
    """
    2 Types of input dictionary are available in Aligner

    * THEMATIC = "thematic"
    * REFERENCE = "reference"
    """

    THEMATIC = "thematic"
    REFERENCE = "reference"


class DiffMetric(str, Enum):
    """
    Metrics to quantify the change between thematic and reference data.

    Attributes
    ----------
    AREA_CHANGE : str
        Absolute change in area (m²).
    AREA_PERCENTAGE_CHANGE : str
        Relative change in area (%).
    SYMMETRICAL_AREA_CHANGE : str
        The area of the symmetric difference (XOR).
    SYMMETRICAL_AREA_PERCENTAGE_CHANGE : str
        Symmetric difference as a percentage of original area.
    LENGTH_CHANGE : str
        Change in perimeter/length.
    LENGTH_PERCENTAGE_CHANGE : str
        Relative change in perimeter.
    LENGTH_REMOVED : str
        The length of segments removed during processing.
    TOTAL_DISTANCE : str
        Sum of displacement of all vertices.
    REFERENCE_USAGE : str
        The extent of reference boundaries utilized (m or m²).
    """

    AREA_CHANGE = "AREA_CHANGE"
    AREA_PERCENTAGE_CHANGE = "AREA_PERCENTAGE_CHANGE"
    SYMMETRICAL_AREA_CHANGE = "SYMMETRICAL_AREA_CHANGE"
    SYMMETRICAL_AREA_PERCENTAGE_CHANGE = "SYMMETRICAL_AREA_PERCENTAGE_CHANGE"
    LENGTH_CHANGE = "LENGTH_CHANGE"
    LENGTH_PERCENTAGE_CHANGE = "LENGTH_PERCENTAGE_CHANGE"
    LENGTH_REMOVED = "LENGTH_REMOVED"
    TOTAL_DISTANCE = "TOTAL_DISTANCE"
    REFERENCE_USAGE = "REFERENCE_USAGE"


class ProcessRemark(str, Enum):
    """
    Status remarks added to processed features for auditing and debugging.

    Attributes
    ----------
    RESULT_UNCHANGED : str
        The alignment did not change the geometry (within tolerance).
    INPUT_CIRCLE : str
        A circle was detected; input returned to prevent geometric distortion.
    RESULT_EMPTY_RETURNED : str
        The processing resulted in an empty geometry.
    CHANGED_GEOMETRYTYPE_EMPTY_RETURNED : str
        Processing changed the geometry type (e.g. Polygon to Line), which is treated as invalid.
    CHANGED_AMOUNT_GEOMETRIES : str
        The number of parts in a Multi-geometry has changed.
    NO_PREDICTION_ORIGINAL_RETURNED : str
        No suitable prediction found; falling back to original.
    MULTIPLE_PREDICTIONS_ORIGINAL_RETURNED : str
        Ambiguous result due to multiple predictions.
    NOT_EVALUATED_ORIGINAL_RETURNED : str
        Evaluation criteria not met.
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
    Classification of the alignment quality and reliability.

    These values help in deciding whether a change can be accepted automatically
    or requires manual review.

    Attributes
    ----------
    EQUALITY_BY_ID_AND_FULL_REFERENCE : str
        Match confirmed by both unique identifier and complete reference data.
    EQUALITY_BY_ID : str
        Match confirmed via unique identifier only.
    EQUALITY_BY_FULL_REFERENCE : str
        Match confirmed via complete reference data comparison.
    PREDICTION_UNIQUE : str
        A single high-confidence prediction was found.
    PREDICTION_UNIQUE_AND_FULL_REFERENCE : str
        A unique prediction that also matches the full reference.
    TO_CHECK_PREDICTION_FULL : str
        Requires review; prediction exists with full reference but lacks certainty.
    TO_CHECK_PREDICTION_MULTI : str
        Requires review; multiple conflicting predictions were found.
    TO_CHECK_PREDICTION_MULTI_FULL : str
        Requires review; multiple predictions found alongside full reference data.
    TO_CHECK_ORIGINAL : str
        Requires review against the original source data.
    TO_CHECK_NO_PREDICTION : str
        Requires review because no valid prediction could be generated.
    NOT_EVALUATED : str
        The alignment has not yet been processed or evaluated.
    NO_CHANGE : str
        Evaluation complete; no changes were detected or required.
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
    Unique identifiers for the available alignment algorithms.

    References correspond to the internal Aligner documentation and methodology.
    """

    DIEUSSAERT = "2024:dieussaert2024a"
    SNAP = "2024:snap2024a"
    NETWORK = "2024:network2024a"
    ALIGNER = "2024:aligner2024a"
    TOPOLOGY = "2024:topology2024a"


class AlignerResultType(str, Enum):
    """Format of the output dictionary."""

    PREDICTIONS = "predictions"
    EVALUATED_PREDICTIONS = "evaluated_predictions"
    PROCESSRESULTS = "processresults"


class AlignerInputType(str, Enum):
    """Role of the input dataset."""

    THEMATIC = "thematic"
    REFERENCE = "reference"


class SnapStrategy(str, Enum):
    """Geometric priority during snapping."""

    ONLY_VERTICES = "only_vertices"
    PREFER_VERTICES = "prefer_vertices"
    NO_PREFERENCE = "no_preference"


class FullReferenceStrategy(str, Enum):
    """Strategy for handling reference data coverage."""

    ONLY_FULL_REFERENCE = "only_full_reference"
    PREFER_FULL_REFERENCE = "prefer_full_reference"
    NO_FULL_REFERENCE = "no_full_reference"


class PredictionStrategy(str, Enum):
    """Determines which prediction is selected as the primary result."""

    ALL = "all"
    BEST = "best"
    ORIGINAL = "original"
