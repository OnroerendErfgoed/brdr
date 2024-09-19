from enum import Enum
from enum import IntEnum


class OpenbaarDomeinStrategy(IntEnum):
    """
    Determines how the algorithm deals with parts of the geometry that are not on the
    reference layer. (=public domain in the case of plots as a reference layer).
    Different strategies have been developed:

    *   -1 (EXCLUDE): Completely exclude everything that is not on the reference layer
    *   0 (AS_IS): All parts that are not covered by the reference layer are added to
        the resulting geometry AS IS
    *   1 (SNAP_SINGLE_SIDE): Everything that falls within the relevant distance over
        the plot boundary is snapped to the plot. The outer boundary is not used.
    *   2 (SNAP_ALL_SIDE): Everything that falls within the relevant distance over the
        plot boundary is snapped to the plot. The outer boundary is used.
    *   3 (SNAP_FULL_AREA_SINGLE_SIDE): integrates entire Openbaar Domein within
        buffered original geometry
    *   4 (SNAP_FULL_AREA_ALL_SIDE): integrates all OD within inside-buffered
        original geometry. The first part is a copy of SNAP_ALL_SIDE
    *   5 (SNAP_SINGLE_SIDE_VARIANT_1): implementation variant 1 of strategy
        SNAP_SINGLE_SIDE
    *   6 (SNAP_SINGLE_SIDE_VARIANT_2): implementation variant 2 of strategy
        SNAP_SINGLE_SIDE

    """

    EXCLUDE = -1
    AS_IS = 0
    SNAP_SINGLE_SIDE = 1
    SNAP_ALL_SIDE = 2
    SNAP_FULL_AREA_SINGLE_SIDE = 3
    SNAP_FULL_AREA_ALL_SIDE = 4
    SNAP_SINGLE_SIDE_VARIANT_1 = 5
    SNAP_SINGLE_SIDE_VARIANT_2 = 6


class AlignerResultType(str, Enum):
    """
    2 Types of resulting dictionary are available in Aligner

    * PREDICTIONS = "predictions" (only the predicted versions for specific relevenat distances are returned)
    * PROCESSRESULTS = "processresults" (All versions of resulting geometries for all relevant distances are returned)
    """

    PREDICTIONS = "predictions"
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
    """

    TOTAL_AREA = "total_area"
    TOTAL_PERCENTAGE = "total_percentage"
    CHANGES_AREA = "changes_area"
    CHANGES_PERCENTAGE = "changes_percentage"


class Evaluation(str, Enum):
    """
    Enum to evaluate an automatically updated geometry:

    * EQUALITY_FORMULA_GEOM_1 = "equality_formula_geom_1"
    * EQUALITY_FORMULA_2 = "equality_formula_2"
    * EQUALITY_GEOM_3 = "equality_geom_3"
    * TO_CHECK_4 = "to_check_4"
    * NO_PREDICTION_5 = "no_prediction_5"
    * NO_CHANGE_6 = "no_change_6"
    """

    EQUALITY_FORMULA_GEOM_1 = "equality_formula_geom_1"
    EQUALITY_FORMULA_2 = "equality_formula_2"
    EQUALITY_GEOM_3 = "equality_geom_3"
    TO_CHECK_4 = "to_check_4"
    NO_PREDICTION_5 = "no_prediction_5"
    NO_CHANGE_6 = "no_change_6"
