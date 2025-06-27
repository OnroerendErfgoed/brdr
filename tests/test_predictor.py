import unittest

import numpy as np
from shapely import from_wkt
from shapely.geometry import Polygon

from brdr.aligner import Aligner
from brdr.enums import GRBType
from brdr.grb import (
    GRBActualLoader,
)
from brdr.loader import DictLoader


class TestAligner(unittest.TestCase):
    def setUp(self):
        # Create a sample geometry for testing
        self.sample_aligner = Aligner()
        self.sample_geom = Polygon([(0, 0), (0, 10), (10, 10), (10, 0)])

    def test_predictor(self):
        # Load thematic data & reference data
        thematic_dict = {"theme_id": from_wkt("POLYGON ((0 0, 0 9, 5 10, 10 0, 0 0))")}
        # ADD A REFERENCE POLYGON TO REFERENCE DICTIONARY
        reference_dict = {"ref_id": from_wkt("POLYGON ((0 1, 0 10,8 10,10 1,0 1))")}
        # LOAD THEMATIC DICTIONARY
        self.sample_aligner.load_thematic_data(DictLoader(thematic_dict))
        # LOAD REFERENCE DICTIONARY
        self.sample_aligner.load_reference_data(DictLoader(reference_dict))
        series = np.arange(0, 310, 10, dtype=int) / 100
        # predict which relevant distances are interesting to propose as resulting
        # geometry

        dict_series, dict_predictions, dict_diffs = self.sample_aligner.predictor(
            relevant_distances=series, od_strategy=4, threshold_overlap_percentage=50
        )
        self.assertEqual(len(dict_predictions), len(thematic_dict))

    def test_predictor_double_prediction(self):
        """
        Test if a double prediction is filtered out of the prediction results. (from 5 predictions to 4 predictions)
        This testdata has 2 resulting predictions that are the same (at 0.0 and 6.0), and 6.0 will be removed from dict_predictions
        """
        # Initiate an Aligner
        aligner = Aligner()
        # Load thematic data & reference data
        loader = DictLoader(
            {
                "id1": from_wkt(
                    "MultiPolygon Z (((138430.4033999964594841 194082.86080000177025795 0, 138422.19659999758005142 194080.36510000005364418 0, 138419.01550000160932541 194079.34930000081658363 0, 138412.59849999845027924 194077.14139999821782112 0, 138403.65579999983310699 194074.06430000066757202 0, 138402.19910000264644623 194077.67480000108480453 0, 138401.83420000225305557 194078.57939999923110008 0, 138400.89329999685287476 194080.91140000149607658 0, 138400.31650000065565109 194080.67880000174045563 0, 138399.27300000190734863 194083.37680000066757202 0, 138405.93310000002384186 194085.95410000160336494 0, 138413.51049999892711639 194088.80620000138878822 0, 138427.25680000334978104 194094.29969999939203262 0, 138430.4033999964594841 194082.86080000177025795 0)))"
                )
            }
        )
        aligner.load_thematic_data(loader)
        loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
        aligner.load_reference_data(loader)

        series = np.arange(0, 810, 10, dtype=int) / 100
        # predict which relevant distances are interesting to propose as resulting geometry
        dict_series, dict_predictions, diffs = aligner.predictor(
            relevant_distances=series, od_strategy=4, threshold_overlap_percentage=50
        )
        assert len(dict_predictions["id1"]) >= 1

    def test_predictor_no_prediction(self):
        """
        Test if no prediction is returned when there are no stable geometries (=no zerostreaks) in the range of relevent_distances
        This testdata has no prediction between 0 and 1 meter, constantly increasing diff (first stable geometry at relevant_distances>2)
        """
        # Initiate an Aligner
        aligner = Aligner(max_workers=-1)
        # Load thematic data & reference data
        wkt = "Polygon ((174125.66583829143201001 179320.36133541050367057, 174125.66583829294540919 179320.36133543887990527, 174123.56114244047785178 179320.47114562886417843, 174123.60579274909105152 179320.9312379399780184, 174123.65622600910137407 179321.4509197928418871, 174130.00431616476271302 179386.86384879014804028, 174131.20465914410306141 179386.55606853921199217, 174131.47325675669708289 179386.48719735653139651, 174131.4777007331722416 179386.48605787538690493, 174131.48166125148418359 179386.48504236005828716, 174131.56964346842141822 179386.46248281505540945, 174145.93543933291221038 179382.77894541397108696, 174153.73543597472598776 179380.77894627506611869, 174155.22223844396648929 179380.39771487266989425, 174155.31905292189912871 179380.37289064755896106, 174158.08159073328715749 179379.66454761903150938, 174162.78046076380996965 179378.45970914987265132, 174162.77348954943590797 179378.41265345277497545, 174162.76924267943832092 179378.38398708027671091, 174162.76793666195590049 179378.37517146224854514, 174158.29825295542832464 179348.20480644324561581, 174158.06896705977851525 179348.24847994712763466, 174157.82196880917763337 179348.29552723292727023, 174152.4153556079545524 179349.32535831883433275, 174151.38357602406176738 179340.0393457512545865, 174150.45438929076772183 179331.6766651498619467, 174155.32298101272317581 179329.59012584027368575, 174155.97540604253299534 179329.31051516399020329, 174156.06587580699124373 179329.27174254972487688, 174156.3372870393213816 179329.15542325720889494, 174155.77701106332824565 179318.79031770044821315, 174154.96588003114447929 179318.83263741704286076, 174154.95658402171102352 179318.83312241756357253, 174151.97826524809352122 179318.98851313433260657, 174146.5353382924804464 179319.27249193197349086, 174125.66583829143201001 179320.36133541050367057))"

        loader = DictLoader({"id1": from_wkt(wkt)})
        aligner.load_thematic_data(loader)
        loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
        aligner.load_reference_data(loader)

        series = np.arange(0, 110, 10, dtype=int) / 100
        # predict which relevant distances are interesting to propose as resulting geometry
        dict_series, dict_predictions, diffs = aligner.predictor(
            relevant_distances=series, od_strategy=4, threshold_overlap_percentage=50
        )
        assert len(dict_predictions["id1"]) <= 1

    def test_predictor_point(self):
        # Load thematic data & reference data
        thematic_dict = {"theme_id": from_wkt("POINT (0 0)")}
        # ADD A REFERENCE POLYGON TO REFERENCE DICTIONARY
        reference_dict = {"ref_id": from_wkt("POLYGON ((0 1, 0 10,8 10,10 1,0 1))")}
        # LOAD THEMATIC DICTIONARY
        self.sample_aligner.load_thematic_data(DictLoader(thematic_dict))
        # LOAD REFERENCE DICTIONARY
        self.sample_aligner.load_reference_data(DictLoader(reference_dict))
        series = np.arange(0, 310, 10, dtype=int) / 100
        # predict which relevant distances are interesting to propose as resulting
        # geometry

        dict_series, dict_predictions, dict_diffs = self.sample_aligner.predictor(
            relevant_distances=series,
        )
        self.assertEqual(len(dict_predictions), len(thematic_dict))
        assert dict_predictions["theme_id"][0.0]["brdr_prediction_count"] >= 1

    def test_predictor_line(self):
        # Load thematic data & reference data
        thematic_dict = {
            "theme_id": from_wkt(
                "LINESTRING (171741.11190000033820979 171839.01070547936251387, 171751.68948904142598622 171847.23771917796693742, 171762.26707808251376264 171855.69979041084297933, 171772.72713835645117797 171862.9865739724773448, 171783.53978493178146891 171870.8610013697471004, 171794.94007534274714999 171877.79519862998859026, 171801.87427260301774368 171884.2592808217741549)"
            )
        }
        self.sample_aligner.load_thematic_data(DictLoader(thematic_dict))
        # LOAD REFERENCE DICTIONARY
        self.sample_aligner.load_reference_data(
            GRBActualLoader(
                aligner=self.sample_aligner, grb_type=GRBType.ADP, partition=1000
            )
        )
        series = np.arange(0, 310, 10, dtype=int) / 100
        # predict which relevant distances are interesting to propose as resulting
        # geometry
        dict_series, dict_predictions, dict_diffs = self.sample_aligner.predictor(
            relevant_distances=series,
        )
        self.assertEqual(len(dict_predictions), len(thematic_dict))
        assert dict_predictions["theme_id"][0.0]["brdr_prediction_count"] >= 1

    def test_predictor_poly_to_point(self):
        # Load thematic data & reference data
        thematic_dict = {
            "theme_id": from_wkt(
                "MULTIPOLYGON (((171807.01653832942247391 171811.64178536087274551, 171806.50108232349157333 171811.20914536342024803, 171801.84252232313156128 171816.64696936681866646, 171798.38581831753253937 171820.68191336840391159, 171767.47471430152654648 171856.76376939192414284, 171767.35881029814481735 171856.89906539395451546, 171771.46146630495786667 171859.87231339514255524, 171777.14249030500650406 171863.9894334003329277, 171779.40213830769062042 171865.62706540152430534, 171783.0694023072719574 171860.27736939489841461, 171789.40066631138324738 171851.0414653904736042, 171790.78703431785106659 171849.01906538754701614, 171812.97858633100986481 171816.6465853676199913, 171809.78524232655763626 171813.96594536304473877, 171807.01653832942247391 171811.64178536087274551)))"
            )
        }
        self.sample_aligner.load_thematic_data(DictLoader(thematic_dict))
        # LOAD REFERENCE DICTIONARY
        reference_dict = {
            "ref_id_1": from_wkt(
                "POINT (171766.12778689494007267 171857.04121096828021109)"
            ),
            "ref_id_2": from_wkt(
                "POINT (171777.13617191879893653 171865.19089055553195067)"
            ),
        }
        # LOAD THEMATIC DICTIONARY
        self.sample_aligner.load_thematic_data(DictLoader(thematic_dict))
        # LOAD REFERENCE DICTIONARY
        self.sample_aligner.load_reference_data(DictLoader(reference_dict))

        series = np.arange(0, 310, 20, dtype=int) / 100
        # predict which relevant distances are interesting to propose as resulting
        # geometry
        dict_series, dict_predictions, dict_diffs = self.sample_aligner.predictor(
            relevant_distances=series,
        )
        self.assertEqual(len(dict_predictions), len(thematic_dict))
        assert dict_predictions["theme_id"][0.0]["brdr_prediction_count"] >= 1
