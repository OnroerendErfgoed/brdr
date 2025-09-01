import unittest

from brdr.aligner import Aligner
from brdr.geometry_utils import geom_from_wkt
from brdr.loader import DictLoader
from brdr.osm import OSMLoader


class TestOsm(unittest.TestCase):

    def test_osmloader(self):
        aligner = Aligner()
        thematic_dict = {
            "my_building_id": geom_from_wkt(
                "POLYGON ((172355.21254837638116442 171905.62395133438985795, 172365.70286233740625903 171900.59734256137744524, 172359.91133483810699545 171887.48445011011790484, 172357.17948224407155067 171887.70299831763259135, 172355.86819299895432778 171884.42477520479587838, 172346.90771649059024639 171888.14009473266196437, 172355.21254837638116442 171905.62395133438985795))"
            )
        }
        loader = DictLoader(data_dict=thematic_dict)
        aligner.load_thematic_data(loader)
        # Use a OSMLoader for the reference data
        loader = OSMLoader(osm_tags={"building": True}, aligner=aligner)
        aligner.load_reference_data(loader)
        assert len(aligner.dict_reference.keys()) > 0
