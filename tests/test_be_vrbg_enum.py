from brdr.be.enums import VRBGType


def test_vrbg_enum_values():
    assert VRBGType.REFGEM.value == "Refgem"
    assert VRBGType.REFARR.value == "Refarr"
    assert VRBGType.REFPRV.value == "Refprv"
