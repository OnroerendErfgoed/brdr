import unittest
from datetime import date, timedelta

from shapely import Polygon, from_wkt

from brdr.aligner import Aligner
from brdr.constants import EVALUATION_FIELD_NAME
from brdr.enums import GRBType, Evaluation
from brdr.grb import (
    get_last_version_date,
    is_grb_changed,
    get_geoms_affected_by_grb_change,
    GRBSpecificDateParcelLoader,
    update_to_actual_grb,
)
from brdr.loader import DictLoader


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
        aligner = Aligner()
        aligner.load_thematic_data(DictLoader(thematic_dict))
        dict_affected, dict_unchanged = get_geoms_affected_by_grb_change(
            aligner=aligner,
            grb_type=GRBType.ADP,
            date_start=date.today() - timedelta(days=30),
            date_end=date.today(),
            one_by_one=False,
            border_distance=0,
        )
        assert len(dict_affected.keys()) > 0

        dict_affected, dict_unchanged = get_geoms_affected_by_grb_change(
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
                "MultiPolygon (((174184.09476602054201066 171899.68933439542888664, "
                "174400.56834639035514556 171832.959863749332726, "
                "174388.65236948925303295 171770.99678386366576888, "
                "174182.10876987033407204 171836.13745758961886168, "
                "174184.88916448061354458 171873.07698598300339654, "
                "174184.09476602054201066 171899.68933439542888664)))"
            )
        }
        aligner = Aligner()
        aligner.load_thematic_data(DictLoader(thematic_dict))
        dict_affected, dict_unchanged = get_geoms_affected_by_grb_change(
            aligner=aligner,
            grb_type=GRBType.ADP,
            date_start=date.today() - timedelta(days=1),
            date_end=date.today(),
            one_by_one=True,
        )
        assert len(dict_affected.keys()) == 0

        dict_affected, dict_unchanged = get_geoms_affected_by_grb_change(
            aligner=aligner,
            grb_type=GRBType.ADP,
            date_start=date.today() - timedelta(days=1000),
            date_end=date.today(),
            one_by_one=True,
        )
        assert len(dict_affected.keys()) == 0
        thematic_dict2 = {
            "theme_id_2": from_wkt(
                "MultiPolygon (((174180.20077791667426936 171966.14649116666987538, "
                "174415.60530965600628406 171940.9636807945498731, "
                "174388.65236948925303295 171770.99678386366576888, "
                "174182.10876987033407204 171836.13745758961886168, "
                "174184.88916448061354458 171873.07698598300339654, "
                "174180.20077791667426936 171966.14649116666987538)))"
            )
        }
        aligner2 = Aligner()
        aligner2.load_thematic_data(DictLoader(thematic_dict2))
        dict_affected, dict_unchanged = get_geoms_affected_by_grb_change(
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
                "MultiPolygon (((174180.20077791667426936 171966.14649116666987538, "
                "174415.60530965600628406 171940.9636807945498731, "
                "174388.65236948925303295 171770.99678386366576888, "
                "174182.10876987033407204 171836.13745758961886168, "
                "174184.88916448061354458 171873.07698598300339654, "
                "174180.20077791667426936 171966.14649116666987538)))"
            )
        }
        aligner = Aligner()
        aligner.load_thematic_data(DictLoader(thematic_dict))
        dict_affected, dict_unchanged = get_geoms_affected_by_grb_change(
            aligner=aligner,
            grb_type=GRBType.ADP,
            date_start=date.today() - timedelta(days=1),
            date_end=date.today(),
            one_by_one=False,
        )
        assert len(dict_affected.keys()) == 0

        dict_affected, dict_unchanged = get_geoms_affected_by_grb_change(
            aligner=aligner,
            grb_type=GRBType.ADP,
            date_start=date.today() - timedelta(days=1000),
            date_end=date.today(),
            one_by_one=False,
        )
        assert len(dict_affected.keys()) > 0

    def test_grbspecificdateparcelloader(self):
        thematic_dict = {
            "theme_id_1": from_wkt(
                "Polygon ((172283.76869662097305991 174272.85233648214489222, 172276.89871930953813717 174278.68436246179044247, 172274.71383684969623573 174280.57171753142029047, 172274.63047763772192411 174280.64478165470063686, 172272.45265833073062822 174282.52660570573061705, 172269.33533191855531186 174285.22093996312469244, 172265.55258252174826339 174288.49089696351438761, 172258.77032718938426115 174294.22654021997004747, 172258.63259260458289646 174294.342757155187428, 172254.93673790179309435 174288.79932878911495209, 172248.71360730109154247 174279.61860501393675804, 172248.96566232520854101 174279.43056782521307468, 172255.25363882273086347 174274.73737183399498463, 172257.08298882702365518 174273.37133203260600567, 172259.32325354730710387 174271.69890458136796951, 172261.65807284769834951 174269.9690355472266674, 172266.35596220899606124 174266.4871726930141449, 172273.34350050613284111 174261.30863015633076429, 172289.60360219911672175 174249.35944479051977396, 172293.30328181147342548 174246.59864199347794056, 172297.34760522318538278 174253.10583685990422964, 172289.53060952731175348 174259.6846851697191596, 172292.86485871637705714 174265.19099397677928209, 172283.76869662097305991 174272.85233648214489222))"
            )
        }
        aligner = Aligner()
        loader = DictLoader(thematic_dict)
        aligner.load_thematic_data(loader)
        loader = GRBSpecificDateParcelLoader(date="2023-01-03", aligner=aligner)
        aligner.load_reference_data(loader)
        assert len(aligner.dict_reference.keys()) == 53

        loader = GRBSpecificDateParcelLoader(date="2023-08-03", aligner=aligner)
        aligner.load_reference_data(loader)
        assert len(aligner.dict_reference.keys()) == 52

    def test_update_to_actual_grb(self):
        # Create a featurecollection (aligned on 2022), to use for the 'update_to_actual_grb'
        name_thematic_id = "theme_identifier"
        featurecollection_base_result = {
            "crs": {"properties": {"name": "EPSG:31370"}, "type": "name"},
            "features": [
                {
                    "geometry": {
                        "coordinates": [
                            [
                                [138541.4173, 194007.7292],
                                [138539.3263, 193994.1382],
                                [138529.3663, 193995.5664],
                                [138522.0997, 193996.6084],
                                [138514.9844, 193997.6287],
                                [138505.8261, 193996.615],
                                [138498.8406, 193996.4314],
                                [138492.9442, 193996.2895],
                                [138491.2246, 193996.2481],
                                [138491.4111, 194004.8147],
                                [138514.3685, 194005.1297],
                                [138520.2585, 194004.5753],
                                [138520.3946, 194005.5833],
                                [138520.5426, 194009.732],
                                [138541.4173, 194007.7292],
                            ]
                        ],
                        "type": "Polygon",
                    },
                    "properties": {
                        "area": 503.67736346047064,
                        "brdr_formula": '{"alignment_date": "2024-09-19", "brdr_version": "0.2.1", "reference_source": {"source": "Adpf", "version_date": "2022-01-01"}, "full": true, "reference_features": {"12034A0181/00K000": {"full": true, "area": 503.68, "percentage": 100, "version_date": "2019-08-30"}}, "reference_od": null, "last_version_date": "2019-08-30"}',
                        "nr_calculations": 1,
                        "perimeter": 125.74541473322422,
                        "relevant_distance": 2,
                        "shape_index": 0.24965468741597102,
                        "theme_identifier": "206285",
                    },
                    "type": "Feature",
                }
            ],
            "type": "FeatureCollection",
        }

        # Update Featurecollection to actual version
        featurecollection = update_to_actual_grb(
            featurecollection_base_result, name_thematic_id
        )
        # Print results
        for feature in featurecollection["result"]["features"]:
            assert isinstance(feature["properties"][EVALUATION_FIELD_NAME], Evaluation)
