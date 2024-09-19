import unittest

from brdr.aligner import Aligner
from brdr.oe import OnroerendErfgoedLoader, OEType


class TestOE(unittest.TestCase):
    def test_onroerenderfgoedloader_by_aanduidid(self):
        loader = OnroerendErfgoedLoader(objectids=[120288, 10275], oetype=OEType.AO)
        aligner = Aligner()
        aligner.load_thematic_data(loader)
        assert len(aligner.dict_thematic.keys()) == 2

    def test_onroerenderfgoedloader_by_erfgoedid(self):
        loader = OnroerendErfgoedLoader(objectids=[42549], oetype=OEType.EO)
        aligner = Aligner()
        aligner.load_thematic_data(loader)
        assert len(aligner.dict_thematic.keys()) == 1

    def test_onroerenderfgoedloader_by_bbox(self):
        loader = OnroerendErfgoedLoader(
            bbox=[172000, 172000, 174000, 174000], oetype=OEType.EO
        )
        aligner = Aligner()
        aligner.load_thematic_data(loader)
        assert len(aligner.dict_thematic.keys()) > 0

    def test_onroerenderfgoedloader_by_bbox_and_objectid(self):

        with self.assertRaises(Exception) as context:
            loader = OnroerendErfgoedLoader(
                objectids=[42549],
                bbox=[172000, 172000, 174000, 174000],
                oetype=OEType.EO,
            )

        with self.assertRaises(Exception) as context:
            loader = OnroerendErfgoedLoader(objectids=None, bbox=None, oetype=OEType.EO)
