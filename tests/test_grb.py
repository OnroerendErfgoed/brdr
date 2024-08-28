import unittest
from datetime import date, timedelta

import numpy as np
from shapely import Polygon, from_wkt

from brdr.aligner import Aligner
from brdr.enums import GRBType
from brdr.geometry_utils import get_bbox
from brdr.grb import (
    get_last_version_date,
    is_grb_changed,
    get_geoms_affected_by_grb_change, evaluate, get_collection_grb_fiscal_parcels,
)
from brdr.loader import DictLoader, GeoJsonLoader
from brdr.loader import GRBActualLoader
from brdr.utils import (get_series_geojson_dict,
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
        aligner=Aligner()
        aligner.load_thematic_data(DictLoader(thematic_dict))
        dict_affected = get_geoms_affected_by_grb_change(
            aligner=aligner,
            grb_type=GRBType.ADP,
            date_start=date.today() - timedelta(days=30),
            date_end=date.today(),
            one_by_one=False,
            border_distance=0,
        )
        assert len(dict_affected.keys()) > 0

        dict_affected = get_geoms_affected_by_grb_change(
            aligner=aligner,
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
        aligner=Aligner()
        aligner.load_thematic_data(DictLoader(thematic_dict))
        dict_affected = get_geoms_affected_by_grb_change(
            aligner=aligner,
            grb_type=GRBType.ADP,
            date_start=date.today() - timedelta(days=1),
            date_end=date.today(),
            one_by_one=True,
        )
        assert len(dict_affected.keys()) == 0

        dict_affected = get_geoms_affected_by_grb_change(
            aligner=aligner,
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
        aligner2=Aligner()
        aligner2.load_thematic_data(DictLoader(thematic_dict2))
        dict_affected = get_geoms_affected_by_grb_change(
            aligner=aligner2,
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
        aligner=Aligner()
        aligner.load_thematic_data(DictLoader(thematic_dict))
        dict_affected = get_geoms_affected_by_grb_change(
            aligner=aligner,
            grb_type=GRBType.ADP,
            date_start=date.today() - timedelta(days=1),
            date_end=date.today(),
            one_by_one=False,
        )
        assert len(dict_affected.keys()) == 0

        dict_affected = get_geoms_affected_by_grb_change(
            aligner=aligner,
            grb_type=GRBType.ADP,
            date_start=date.today() - timedelta(days=1000),
            date_end=date.today(),
            one_by_one=False,
        )
        assert len(dict_affected.keys()) > 0

    def test_evaluate(self):
        thematic_dict = {
            "theme_id_1": from_wkt(
                "MultiPolygon (((174180.20077791667426936 171966.14649116666987538, 174415.60530965600628406 171940.9636807945498731, 174388.65236948925303295 171770.99678386366576888, 174182.10876987033407204 171836.13745758961886168, 174184.88916448061354458 171873.07698598300339654, 174180.20077791667426936 171966.14649116666987538)))"
            )
        }
        bbox = get_bbox(thematic_dict["theme_id_1"])
        base_aligner = Aligner()
        base_aligner.load_thematic_data(DictLoader(thematic_dict))
        base_year = "2022"
        collection_fiscal_parcels = get_collection_grb_fiscal_parcels(base_year, bbox=bbox)
        base_aligner.load_reference_data(GeoJsonLoader(collection_fiscal_parcels, "CAPAKEY"))
        base_process_result = base_aligner.process_dict_thematic(relevant_distance=1)
        thematic_dict_formula = {}
        thematic_dict_result = {}
        for key in base_process_result:
            thematic_dict_result[key] = base_process_result[key]["result"]
            thematic_dict_formula[key] = base_aligner.get_formula(thematic_dict_result[key])
        aligner_result=Aligner()
        aligner_result.load_thematic_data(DictLoader(thematic_dict_result))
        dict_affected = get_geoms_affected_by_grb_change(
            aligner=aligner_result,
            grb_type=GRBType.ADP,
            date_start=date(2022, 1, 1),
            date_end=date.today(),
            one_by_one=False,
        )

        series = np.arange(0, 200, 10, dtype=int) / 100

        actual_aligner = Aligner()
        loader = DictLoader(dict_affected)
        actual_aligner.load_thematic_data(loader)
        loader = GRBActualLoader(grb_type=GRBType.ADP, partition=0, aligner=actual_aligner)
        actual_aligner.load_reference_data(loader)

        dict_evaluated, prop_dictionary = evaluate(actual_aligner, thematic_dict_formula, series, )
        fc = get_series_geojson_dict(
            dict_evaluated,
            crs=actual_aligner.CRS,
            id_field=actual_aligner.name_thematic_id,
            series_prop_dict=prop_dictionary,
        )

        print(fc)


