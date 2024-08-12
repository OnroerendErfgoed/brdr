import unittest
from datetime import date, datetime, timedelta

from shapely import Polygon, from_wkt

from brdr.aligner import Aligner
from brdr.enums import GRBType
from brdr.grb import (
    get_last_version_date,
    is_grb_changed,
    get_geoms_affected_by_grb_change,
)
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

    def test_is_grb_changed_outerborder(self):
        geom = Polygon(
            [(170000, 170000), (170000, 172000), (172000, 172000), (172000, 170000)]
        )
        out = is_grb_changed(
            geom,
            border_distance=0,
            grb_type=GRBType.ADP,
            date_start=date(2024, 7, 1),
        )
        self.assertTrue(out)
        out = is_grb_changed(
            geom,
            border_distance=10,
            grb_type=GRBType.ADP,
            date_start=date(2024, 7, 1),
        )
        self.assertFalse(out)

    def test_get_geoms_affected_by_grb_change_outerborder(self):
        thematic_dict = {
            "theme_id_1": Polygon(
                [(170000, 170000), (170000, 172000), (172000, 172000), (172000, 170000)]
            )
        }

        dict_affected = get_geoms_affected_by_grb_change(
            thematic_dict,
            grb_type=GRBType.ADP,
            date_start=date.today() - timedelta(days=30),
            date_end=date.today(),
            one_by_one=False,
            border_distance=0,
        )
        assert len(dict_affected.keys()) > 0

        dict_affected = get_geoms_affected_by_grb_change(
            thematic_dict,
            grb_type=GRBType.ADP,
            date_start=date.today() - timedelta(days=30),
            date_end=date.today(),
            one_by_one=False,
            border_distance=10,
        )
        assert len(dict_affected.keys()) == 0

    def test_get_geoms_affected_by_grb_change(self):
        thematic_dict = {
            "theme_id_1": from_wkt(
                "MultiPolygon (((174184.09476602054201066 171899.68933439542888664, 174400.56834639035514556 171832.959863749332726, 174388.65236948925303295 171770.99678386366576888, 174182.10876987033407204 171836.13745758961886168, 174184.88916448061354458 171873.07698598300339654, 174184.09476602054201066 171899.68933439542888664)))"
            )
        }
        dict_affected = get_geoms_affected_by_grb_change(
            thematic_dict,
            grb_type=GRBType.ADP,
            date_start=date.today() - timedelta(days=1),
            date_end=date.today(),
            one_by_one=True,
        )
        assert len(dict_affected.keys()) == 0

        dict_affected = get_geoms_affected_by_grb_change(
            thematic_dict,
            grb_type=GRBType.ADP,
            date_start=date.today() - timedelta(days=1000),
            date_end=date.today(),
            one_by_one=True,
        )
        assert len(dict_affected.keys()) == 0
        thematic_dict2 = {
            "theme_id_2": from_wkt(
                "MultiPolygon (((174180.20077791667426936 171966.14649116666987538, 174415.60530965600628406 171940.9636807945498731, 174388.65236948925303295 171770.99678386366576888, 174182.10876987033407204 171836.13745758961886168, 174184.88916448061354458 171873.07698598300339654, 174180.20077791667426936 171966.14649116666987538)))"
            )
        }
        dict_affected = get_geoms_affected_by_grb_change(
            thematic_dict2,
            grb_type=GRBType.ADP,
            date_start=date.today() - timedelta(days=1000),
            date_end=date.today(),
            one_by_one=True,
        )
        assert len(dict_affected.keys()) > 0

    def test_get_geoms_affected_by_grb_change_bulk(self):
        thematic_dict = {
            "theme_id_1": from_wkt(
                "MultiPolygon (((174180.20077791667426936 171966.14649116666987538, 174415.60530965600628406 171940.9636807945498731, 174388.65236948925303295 171770.99678386366576888, 174182.10876987033407204 171836.13745758961886168, 174184.88916448061354458 171873.07698598300339654, 174180.20077791667426936 171966.14649116666987538)))"
            )
        }
        dict_affected = get_geoms_affected_by_grb_change(
            thematic_dict,
            grb_type=GRBType.ADP,
            date_start=date.today() - timedelta(days=1),
            date_end=date.today(),
            one_by_one=False,
        )
        assert len(dict_affected.keys()) == 0

        dict_affected = get_geoms_affected_by_grb_change(
            thematic_dict,
            grb_type=GRBType.ADP,
            date_start=date.today() - timedelta(days=1000),
            date_end=date.today(),
            one_by_one=False,
        )
        assert len(dict_affected.keys()) > 0
