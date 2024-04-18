from enum import IntEnum
from enum import StrEnum


class OpenbaarDomeinStrategy(IntEnum):
    """
    Determines how the algorithm deals with parts of the geometry that are not on the
    reference layer. (=public domain in the case of plots as a reference layer).
    Different strategies have been developed:

    *   -1 (EXCLUDE): Completely exclude everything that is not on the reference layer
    *   0 (AS_IS): All parts that are not covered by the reference layer are added to the
        resulting geometry AS IS
    *   1 (SNAP_ONE_SIDED): Everything that falls within the relevant distance over the
        plot boundary is snapped to the plot. The outer boundary is not used.
    *   2 (SNAP_BOTH_SIDED): Everything that falls within the relevant distance over the
        plot boundary is snapped to the plot. The outer boundary is used.
    *   3 (SNAP_BIG_AREA): integrates entire Openbaar Domein within buffered original
        geometry
    *   4 (SNAP_FULL_AREA): integrates all OD within inside-buffered original geometry.
        The first part is a copy of SNAP_BOTH_SIDED
    *   5 (SNAP_FULL_AREA_BOTH_SIDED): integrates all OD within outside-buffered original

    """

    EXCLUDE = -1
    AS_IS = 0
    SNAP_ONE_SIDED = 1
    SNAP_BOTH_SIDED = 2
    SNAP_FULL_AREA_ONE_SIDED = 3
    SNAP_FULL_AREA_BOTH_SIDED = 4


class GRBType(StrEnum):
    """
    Determines which GRB feature collection is used. Different types are available:

    * ADP: administrative plots
    * GBG: buildings on the ground
    * KNW: artworks
    """

    ADP = "adp"
    GBG = "gbg"
    KNW = "knw"
