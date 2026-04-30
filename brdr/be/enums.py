from enum import Enum


class VRBGType(str, Enum):
    """
    Supported VRBG administrative boundary collections.
    """

    GRENSINFO = "Grensinfo"
    REFARR = "Refarr"
    REFGEM = "Refgem"
    REFGEW = "Refgew"
    REFPRV = "Refprv"

