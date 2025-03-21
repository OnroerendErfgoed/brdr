from enum import Enum
from enum import IntEnum


class OpenDomainStrategy(IntEnum):
    """
    Determines how the algorithm deals with parts of the geometry that are not on the
    reference layer. (=public domain in the case of plots as a reference layer).
    Different strategies have been developed:

    *   -1 (EXCLUDE): Completely exclude everything that is not on the reference layer
    *   0 (AS_IS): All parts that are not covered by the reference layer are added to
        the resulting geometry AS IS
    *   1 (SNAP_INNER_SIDE): Everything that falls within the relevant distance over
        the plot boundary is snapped to the plot. The outer boundary is not used.
    *   2 (SNAP_ALL_SIDE): Everything that falls within the relevant distance over the
        plot boundary is snapped to the plot. The inner and outer boundary is used where possible.
    *   3 (SNAP_PREFER_VERTICES): The part on the OD is 'snapped' to the closest reference-polygons.
        Vertices of the reference-polygons are preferred above edges if they are within the relevant distance
    *   4 (SNAP_NO_PREFERENCE): The part on the OD is 'snapped' to the closest reference-polygons.
        The full edge of the reference-polygons is used. (No preference of reference-vertices.
    *   5 (SNAP_ONLY_VERTICES): The part on the OD is 'snapped' to the vertices of reference-polygons.

    """

    EXCLUDE = -1
    AS_IS = 0
    SNAP_INNER_SIDE = 1
    SNAP_ALL_SIDE = 2
    SNAP_PREFER_VERTICES = 3
    SNAP_NO_PREFERENCE = 4
    SNAP_ONLY_VERTICES = 5


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


class GRBType(str, Enum):
    """
    Determines which GRB feature collection is used. Different types are available:

    * ADP: administrative plots
    * GBG: buildings on the ground
    * KNW: artworks
    """

    ADP = "adp"
    GBG = "gbg"
    KNW = "knw"


class DiffMetric(str, Enum):
    """
    Determines which metric is used to determine the difference between the thematic and
    reference data. Different metrics are available:

    * TOTAL_AREA: the total difference area
    * TOTAL_PERCENTAGE: the percentage of the total difference area
    * CHANGES_AREA: the sum of the negative and positive difference areas
    * CHANGES_PERCENTAGE: the percentage of the changed area
    * TOTAL_LENGTH ="total_length"
    * CHANGES_LENGTH = "changes_length"
    * TOTAL_DISTANCE = "total_distance"
    """

    TOTAL_AREA = "total_area"
    TOTAL_PERCENTAGE = "total_percentage"
    CHANGES_AREA = "changes_area"
    CHANGES_PERCENTAGE = "changes_percentage"
    TOTAL_LENGTH = "total_length"
    CHANGES_LENGTH = "changes_length"
    TOTAL_DISTANCE = "total_distance"


class Evaluation(str, Enum):
    """
    Enum to evaluate an automatically updated geometry:

    EQUALITY_EQUAL_FORMULA_FULL_1 = "equality_equal_formula_full_1"
    EQUALITY_EQUAL_FORMULA_2 = "equality_equal_formula_2"
    EQUALITY_FULL_3 = "equality_full_3"
    PREDICTION_UNIQUE = "prediction_unique"
    PREDICTION_UNIQUE_FULL = "prediction_unique_full"
    TO_CHECK_PREDICTION_FULL = "to_check_prediction_full"
    TO_CHECK_PREDICTION_MULTI = "to_check_prediction_multi"
    TO_CHECK_PREDICTION_MULTI_FULL = "to_check_prediction_multi_full"
    TO_CHECK_ORIGINAL ="to_check_original"
    TO_CHECK_NO_PREDICTION = "to_check_no_prediction"
    NO_CHANGE = "no_change"
    """

    EQUALITY_EQUAL_FORMULA_FULL_1 = "equality_equal_formula_full_1"
    EQUALITY_EQUAL_FORMULA_2 = "equality_equal_formula_2"
    EQUALITY_FULL_3 = "equality_full_3"
    PREDICTION_UNIQUE = "prediction_unique"
    PREDICTION_UNIQUE_FULL = "prediction_unique_full"
    TO_CHECK_PREDICTION_FULL = "to_check_prediction_full"
    TO_CHECK_PREDICTION_MULTI = "to_check_prediction_multi"
    TO_CHECK_PREDICTION_MULTI_FULL = "to_check_prediction_multi_full"
    TO_CHECK_ORIGINAL = "to_check_original"
    TO_CHECK_NO_PREDICTION = "to_check_no_prediction"
    NO_CHANGE = "no_change"


class FullStrategy(str, Enum):
    """
    Enum for full strategy when evaluating predictions:

    ONLY_FULL = "only_full"
    PREFER_FULL = "prefer_full"
    NO_FULL = "no_full"
    """

    ONLY_FULL = "only_full"
    PREFER_FULL = "prefer_full"
    NO_FULL = "no_full"


class SnapStrategy(str, Enum):
    """
    Enum for snapping strategy when snapping a polygon to a reference:

    ONLY_VERTICES = "only_vertices"
    PREFER_VERTICES = "prefer_vertices"
    NO_PREFERENCE = "no_preference"
    """

    ONLY_VERTICES = "only_vertices"
    PREFER_VERTICES = "prefer_vertices"
    NO_PREFERENCE = "no_preference"
