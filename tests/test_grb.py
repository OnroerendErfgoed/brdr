import unittest
from datetime import date

from shapely import Polygon

from brdr.aligner import Aligner
from brdr.enums import GRBType
from brdr.grb import get_last_version_date, is_grb_changed
from brdr.loader import DictLoader
from brdr.loader import GRBActualLoader
from brdr.utils import (
    get_oe_dict_by_ids,
)


class TestGrb(unittest.TestCase):
    def test_get_last_version_date(self):
        # Check if the result of the _buffer_neg_pos gives an equal geometry
        geom = Polygon([(0, 0), (0, 10), (10, 10), (10, 0)])
        out = get_last_version_date(geom)
        self.assertIsNone(out)
        geom = Polygon(
            [(170000, 170000), (170000, 170100), (170100, 170100), (170100, 170000)]
        )
        out = get_last_version_date(geom)
        self.assertIsInstance(out, date)

    def test_is_grb_changed(self):
        # Check if the result of the _buffer_neg_pos gives an equal geometry
        geom = Polygon([(0, 0), (0, 10), (10, 10), (10, 0)])
        out = is_grb_changed(geom)
        self.assertIsNone(out)
        geom = Polygon(
            [(170000, 170000), (170000, 170100), (170100, 170100), (170100, 170000)]
        )
        out = is_grb_changed(geom, grb_type=GRBType.ADP, date_start=date(2024, 7, 16))
        self.assertFalse(out)
        out = is_grb_changed(geom, grb_type=GRBType.ADP, date_start=date(2021, 7, 16))
        self.assertTrue(out)

    def test_get_geoms_affected_by_grb_change(self):
        # TODO test
        assert True
