import unittest
from datetime import date, timedelta

from shapely import Polygon, from_wkt

from brdr.aligner import Aligner
from brdr.enums import GRBType, Evaluation
from brdr.grb import (
    get_last_version_date,
    is_grb_changed,
    get_geoms_affected_by_grb_change,
    GRBSpecificDateParcelLoader, update_to_actual_grb,
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
        aligner = Aligner()
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
        assert len (aligner.dict_reference.keys())==53

        loader = GRBSpecificDateParcelLoader(date="2023-08-03", aligner=aligner)
        aligner.load_reference_data(loader)
        assert len (aligner.dict_reference.keys())==52

    def test_grbspecificdateparcelloader(self):
        #Create a featurecollection (aligned on 2022), to use for the 'update_to_actual_grb'
        name_thematic_id = "theme_identifier"
        featurecollection_base_result = {"crs": {"properties": {"name": "EPSG:31370"}, "type": "name"}, "features": [{"geometry": {"coordinates": [[[174165.099014, 179510.530095], [174165.8317, 179512.9879], [174171.989, 179533.6401], [174176.4529, 179548.8062], [174179.309, 179558.51], [174179.380292, 179558.485703], [174181.1589, 179557.8801], [174187.9589, 179555.5901], [174190.259, 179554.81], [174197.229, 179552.4601], [174199.5291, 179551.6901], [174203.588398, 179550.315901], [174204.019, 179550.1701], [174196.945502, 179518.200008], [174196.899, 179517.9901], [174193.6237, 179503.5462], [174193.0752, 179501.1272], [174193.069, 179501.1002], [174192.963218, 179501.135794], [174183.549015, 179504.310095], [174174.279, 179507.4301], [174167.8091, 179509.6149], [174165.099014, 179510.530095]]], "type": "Polygon"}, "properties": {"area": 1390.3280890476424, "formula": "{\"alignment_date\": \"2024-09-16\", \"brdr_version\": \"0.2.1\", \"reference_source\": {\"source\": \"Adpf\", \"version_date\": \"2022-01-01\"}, \"full\": true, \"reference_features\": {\"24126B0049/00X000\": {\"full\": true, \"area\": 502.91, \"percentage\": 100, \"geometry\": null}, \"24126B0049/00Z000\": {\"full\": true, \"area\": 398.32, \"percentage\": 100, \"geometry\": null}, \"24126B0049/00Y000\": {\"full\": true, \"area\": 489.09, \"percentage\": 100, \"geometry\": null}}, \"reference_od\": null, \"last_version_date\": \"2022-07-29\"}", "perimeter": 155.9132823823815, "relevant_distance": 2, "shape_index": 0.11214135973414749, "theme_identifier": "100"}, "type": "Feature"}, {"geometry": {"coordinates": [[[174149.124298, 179571.446101], [174149.4742, 179571.3366], [174140.7496, 179544.3599], [174140.0649, 179544.0909], [174131.8684, 179521.8687], [174127.3538, 179523.3958], [174125.1598, 179524.1334], [174118.177, 179526.5181], [174117.5579, 179526.7295], [174121.3028, 179537.5797], [174134.5641, 179576.001], [174141.4845, 179573.8361], [174149.124298, 179571.446101]]], "type": "Polygon"}, "properties": {"area": 818.9938386019529, "formula": "{\"alignment_date\": \"2024-09-16\", \"brdr_version\": \"0.2.1\", \"reference_source\": {\"source\": \"Adpf\", \"version_date\": \"2022-01-01\"}, \"full\": true, \"reference_features\": {\"24126B0051/00W000\": {\"full\": true, \"area\": 419.99, \"percentage\": 100, \"geometry\": null}, \"24126B0051/00M002\": {\"full\": true, \"area\": 399.01, \"percentage\": 100, \"geometry\": null}}, \"reference_od\": null, \"last_version_date\": \"2022-07-29\"}", "perimeter": 135.6337116105736, "relevant_distance": 2, "shape_index": 0.16561017338311657, "theme_identifier": "200"}, "type": "Feature"}, {"geometry": {"coordinates": [[[174111.549006, 179153.956005], [174111.5042, 179153.9243], [174110.0614, 179154.1094], [174068.867, 179159.3947], [174068.8661, 179159.4262], [174068.8626, 179159.5573], [174073.7483, 179188.9357], [174120.4387, 179180.3235], [174116.1333, 179157.2025], [174111.549006, 179153.956005]]], "type": "Polygon"}, "properties": {"area": 1344.8114559611831, "formula": "{\"alignment_date\": \"2024-09-16\", \"brdr_version\": \"0.2.1\", \"reference_source\": {\"source\": \"Adpf\", \"version_date\": \"2022-01-01\"}, \"full\": true, \"reference_features\": {\"24126B0031/00N005\": {\"full\": true, \"area\": 1344.81, \"percentage\": 100, \"geometry\": null}}, \"reference_od\": null, \"last_version_date\": \"2019-07-25\"}", "perimeter": 149.60002606562426, "relevant_distance": 2, "shape_index": 0.11124237929598835, "theme_identifier": "300"}, "type": "Feature"}, {"geometry": {"coordinates": [[[174034.4177, 178984.8249], [174030.7603, 178982.3136], [174030.6565, 178982.4711], [174025.7399, 178989.9312], [174018.094404, 178999.593195], [174017.939403, 178999.788996], [174016.3725, 179001.7693], [174018.7192, 179003.659], [174021.115, 179005.7825], [174019.7443, 179007.5141], [174019.7371, 179007.5233], [174015.7101, 179025.628], [174040.6882, 179032.2831], [174037.3194, 178987.071901], [174037.2994, 178986.8036], [174036.3836, 178986.1748], [174034.4177, 178984.8249]]], "type": "Polygon"}, "properties": {"area": 842.1930629252586, "formula": "{\"alignment_date\": \"2024-09-16\", \"brdr_version\": \"0.2.1\", \"reference_source\": {\"source\": \"Adpf\", \"version_date\": \"2022-01-01\"}, \"full\": true, \"reference_features\": {\"24126B0031/00T007\": {\"full\": true, \"area\": 842.19, \"percentage\": 100, \"geometry\": null}}, \"reference_od\": null, \"last_version_date\": \"2019-07-25\"}", "perimeter": 130.58810547796506, "relevant_distance": 2, "shape_index": 0.1550572086457025, "theme_identifier": "400"}, "type": "Feature"}, {"geometry": {"coordinates": [[[173966.389028, 179298.100271], [173965.849202, 179298.315899], [173964.192, 179298.978], [173958.0291, 179301.4402], [173953.8952, 179302.1971], [173948.0517, 179303.2669], [173947.9791, 179303.2803], [173945.8891, 179303.6902], [173911.239422, 179309.581196], [173910.388103, 179309.7266], [173909.9886, 179309.7948], [173905.785701, 179319.638098], [173905.1785, 179321.060399], [173900.5608, 179331.8751], [173900.9241, 179331.8081], [173940.7763, 179325.5153], [173944.092, 179324.9918], [173949.739089, 179324.100202], [173962.1865, 179322.1395], [173966.131594, 179321.518001], [173966.499, 179321.4602], [173970.4676, 179319.931], [173974.1291, 179318.5202], [173972.379009, 179313.840224], [173968.2391, 179302.7402], [173968.229604, 179302.716411], [173966.389028, 179298.100271]]], "type": "Polygon"}, "properties": {"area": 1379.498322959166, "formula": "{\"alignment_date\": \"2024-09-16\", \"brdr_version\": \"0.2.1\", \"reference_source\": {\"source\": \"Adpf\", \"version_date\": \"2022-01-01\"}, \"full\": true, \"reference_features\": {\"24126B0006/00M002\": {\"full\": true, \"area\": 386.55, \"percentage\": 100, \"geometry\": null}, \"24126B0006/00E002\": {\"full\": true, \"area\": 409.64, \"percentage\": 100, \"geometry\": null}, \"24126B0006/00N002\": {\"full\": true, \"area\": 108.75, \"percentage\": 100, \"geometry\": null}, \"24126B0006/00F002\": {\"full\": true, \"area\": 474.56, \"percentage\": 100, \"geometry\": null}}, \"reference_od\": null, \"last_version_date\": \"2021-07-07\"}", "perimeter": 178.54520703582963, "relevant_distance": 2, "shape_index": 0.12942763616619105, "theme_identifier": "500"}, "type": "Feature"}, {"geometry": {"coordinates": [[[174240.361258, 179443.003306], [174240.4272, 179443.1969], [174234.5671, 179445.0969], [174241.3871, 179463.097], [174241.474, 179463.0721], [174244.1019, 179471.6328], [174249.4882, 179469.7988], [174254.26, 179468.16], [174256.144, 179467.513], [174254.5936, 179463.058], [174252.2125, 179456.2165], [174251.3099, 179453.623], [174249.5697, 179448.6229], [174249.0652, 179448.8045], [174248.960701, 179448.502502], [174246.296344, 179440.805126], [174240.361258, 179443.003306]]], "type": "Polygon"}, "properties": {"area": 354.31723731849075, "formula": "{\"alignment_date\": \"2024-09-16\", \"brdr_version\": \"0.2.1\", \"reference_source\": {\"source\": \"Adpf\", \"version_date\": \"2022-01-01\"}, \"full\": false, \"reference_features\": {\"24126B0027/00B002\": {\"full\": false, \"area\": 58.59, \"percentage\": 10.55, \"geometry\": null}, \"24126B0027/00R000\": {\"full\": true, \"area\": 161.68, \"percentage\": 100, \"geometry\": null}, \"24126B0027/00K003\": {\"full\": true, \"area\": 134.05, \"percentage\": 100, \"geometry\": null}}, \"reference_od\": null, \"last_version_date\": \"2019-07-25\"}", "perimeter": 82.71815334885017, "relevant_distance": 2, "shape_index": 0.23345788642649634, "theme_identifier": "600"}, "type": "Feature"}], "type": "FeatureCollection"}

        #Update Featurecollection to actual version
        featurecollection = update_to_actual_grb(featurecollection_base_result,name_thematic_id)
        #Print results
        for feature in featurecollection["result"]["features"]:
            assert isinstance(feature["properties"]["evaluation"],Evaluation)
