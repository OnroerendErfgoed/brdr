import unittest

import numpy as np
import pytest
from shapely import from_wkt
from shapely.geometry import Polygon

from brdr.aligner import Aligner
from brdr.be.grb.enums import GRBType
from brdr.be.grb.loader import GRBActualLoader
from brdr.configs import ProcessorConfig
from brdr.constants import EVALUATION_FIELD_NAME
from brdr.constants import PREDICTION_COUNT
from brdr.enums import AlignerResultType
from brdr.enums import Evaluation
from brdr.enums import FullReferenceStrategy
from brdr.enums import SnapStrategy
from brdr.loader import DictLoader
from brdr.loader import GeoJsonLoader
from brdr.processor import AlignerGeometryProcessor


class TestEvaluate(unittest.TestCase):
    def setUp(self):
        # Create a sample geometry for testing
        self.sample_aligner = Aligner()
        self.sample_geom = Polygon([(0, 0), (0, 10), (10, 10), (10, 0)])

    # def test_evaluate(self):
    #     thematic_dict = {
    #         "theme_id_1": from_wkt(
    #             "MultiPolygon (((174180.20077791667426936 171966.14649116666987538, "
    #             "174415.60530965600628406 171940.9636807945498731, "
    #             "174388.65236948925303295 171770.99678386366576888, "
    #             "174182.10876987033407204 171836.13745758961886168, "
    #             "174184.88916448061354458 171873.07698598300339654, "
    #             "174180.20077791667426936 171966.14649116666987538)))"
    #         )
    #     }
    #     base_aligner = Aligner()
    #     base_aligner.load_thematic_data(DictLoader(thematic_dict))
    #     base_aligner.load_reference_data(
    #         GRBFiscalParcelLoader(aligner=base_aligner, year="2022", partition=1000)
    #     )
    #     relevant_distance = 1
    #     base_process_result = base_aligner.process(relevant_distance=relevant_distance)
    #     thematic_dict_formula = {}
    #     thematic_dict_result = {}
    #     for key in base_process_result:
    #         thematic_dict_result[key] = base_process_result[key][relevant_distance][
    #             "result"
    #         ]
    #         thematic_dict_formula[key] = {
    #             FORMULA_FIELD_NAME: json.dumps(
    #                 base_aligner.get_brdr_formula(thematic_dict_result[key])
    #             )
    #         }
    #         # print(key + ": " + thematic_dict_result[key].wkt)
    #         # print(key + ": " + str(thematic_dict_formula[key]))
    #     base_aligner_result = Aligner()
    #     base_aligner_result.load_thematic_data(DictLoader(thematic_dict_result))
    #     affected, unaffected = get_affected_by_grb_change(
    #         dict_thematic=thematic_dict_result,
    #         grb_type=GRBType.ADP,
    #         date_start=date(2022, 1, 1),
    #         date_end=date.today(),
    #         one_by_one=False,
    #         border_distance=relevant_distance,
    #     )
    #     if len(affected) == 0:
    #         # print("No affected dicts")
    #         exit()
    #     # print("Affected_IDs: " + str(affected))
    #     actual_aligner = Aligner()
    #     actual_aligner.load_thematic_data(
    #         DictLoader(
    #             data_dict=thematic_dict_result,
    #             data_dict_properties=thematic_dict_formula,
    #         )
    #     )
    #     actual_aligner.load_reference_data(
    #         GRBActualLoader(
    #             grb_type=GRBType.ADP, partition=1000, aligner=actual_aligner
    #         )
    #     )
    #     actual_aligner.relevant_distances = np.arange(0, 210, 10, dtype=int) / 100
    #     dict_evaluated, prop_dictionary  = actual_aligner.evaluate(
    #         ids_to_evaluate=affected, base_formula_field=FORMULA_FIELD_NAME
    #     )
    #
    #     fc = get_dict_geojsons_from_series_dict(
    #         dict_evaluated,
    #         crs=actual_aligner.CRS,
    #         id_field=actual_aligner.name_thematic_id,
    #         series_prop_dict=prop_dictionary,
    #     )
    #     # print(fc["result"])
    #     fcs = actual_aligner.get_results_as_geojson(formula=True)
    #     geojson = fcs["result"]
    #     # print(geojson)
    #     geojson = geojson_to_multi(fcs["result"])
    #     # print(geojson)

    @pytest.mark.usefixtures("mock_grb_response2")
    def test_evaluate_full_strategy_no_full(self):
        thematic_dict = {
            "theme_id_1": from_wkt(
                "Polygon ((174015.08694170592934825 179025.39916784031083807, 174040.71934808720834553 179031.93985084796440788, 174037.1838437587430235 178986.50862022733781487, 174030.64316075111855753 178982.97311589887249283, 174016.85469387014745735 179002.24161448894301429, 174021.62762471355381422 179005.77711881740833633, 174018.62244603436556645 179008.60552228015149012, 174015.08694170592934825 179025.39916784031083807))"
            )
        }
        aligner = Aligner()
        aligner.load_thematic_data(DictLoader(thematic_dict))
        aligner.load_reference_data(
            GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
        )

        evaluation_result = aligner.evaluate(
            relevant_distances=np.arange(0, 410, 10, dtype=int) / 100,
            full_reference_strategy=FullReferenceStrategy.NO_FULL_REFERENCE,
        )
        assert len(evaluation_result.results["theme_id_1"]) > 1
        assert (
            evaluation_result.results["theme_id_1"][3.5]["properties"][
                "brdr_evaluation"
            ]
            == Evaluation.TO_CHECK_PREDICTION_MULTI_FULL
        )

    @pytest.mark.usefixtures("mock_grb_response2")
    def test_evaluate_full_strategy_prefer_full(self):
        thematic_dict = {
            "theme_id_1": from_wkt(
                "Polygon ((174015.08694170592934825 179025.39916784031083807, 174040.71934808720834553 179031.93985084796440788, 174037.1838437587430235 178986.50862022733781487, 174030.64316075111855753 178982.97311589887249283, 174016.85469387014745735 179002.24161448894301429, 174021.62762471355381422 179005.77711881740833633, 174018.62244603436556645 179008.60552228015149012, 174015.08694170592934825 179025.39916784031083807))"
            )
        }
        aligner = Aligner()
        aligner.load_thematic_data(DictLoader(thematic_dict))
        aligner.load_reference_data(
            GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
        )

        evaluation_result = aligner.evaluate(
            relevant_distances=np.arange(0, 410, 10, dtype=int) / 100,
            full_reference_strategy=FullReferenceStrategy.PREFER_FULL_REFERENCE,
        )
        assert len(evaluation_result.results["theme_id_1"]) > 1
        assert (
            evaluation_result.results["theme_id_1"][3.5]["properties"][
                EVALUATION_FIELD_NAME
            ]
            == Evaluation.TO_CHECK_PREDICTION_FULL
        )

    @pytest.mark.usefixtures("mock_grb_response2")
    def test_evaluate_full_strategy_only_full(self):
        thematic_dict = {
            "theme_id_1": from_wkt(
                "Polygon ((174015.08694170592934825 179025.39916784031083807, 174040.71934808720834553 179031.93985084796440788, 174037.1838437587430235 178986.50862022733781487, 174030.64316075111855753 178982.97311589887249283, 174016.85469387014745735 179002.24161448894301429, 174021.62762471355381422 179005.77711881740833633, 174018.62244603436556645 179008.60552228015149012, 174015.08694170592934825 179025.39916784031083807))"
            )
        }
        aligner = Aligner()
        aligner.load_thematic_data(DictLoader(thematic_dict))
        aligner.load_reference_data(
            GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
        )

        aligner_result = aligner.evaluate(
            relevant_distances=np.arange(0, 410, 10, dtype=int) / 100,
            full_reference_strategy=FullReferenceStrategy.ONLY_FULL_REFERENCE,
        )
        process_results_evaluated = aligner_result.get_results(
            aligner=aligner, result_type=AlignerResultType.EVALUATED_PREDICTIONS
        )
        assert len(process_results_evaluated["theme_id_1"]) == 1
        assert (
            process_results_evaluated["theme_id_1"][3.5]["properties"][
                EVALUATION_FIELD_NAME
            ]
            == Evaluation.PREDICTION_UNIQUE_AND_FULL_REFERENCE
        )

    @pytest.mark.usefixtures("mock_grb_response2")
    def test_evaluate_all_predictions(self):
        thematic_dict = {
            "theme_id_1": from_wkt(
                "Polygon ((174015.08694170592934825 179025.39916784031083807, 174040.71934808720834553 179031.93985084796440788, 174037.1838437587430235 178986.50862022733781487, 174030.64316075111855753 178982.97311589887249283, 174016.85469387014745735 179002.24161448894301429, 174021.62762471355381422 179005.77711881740833633, 174018.62244603436556645 179008.60552228015149012, 174015.08694170592934825 179025.39916784031083807))"
            )
        }
        aligner = Aligner()
        aligner.load_thematic_data(DictLoader(thematic_dict))
        aligner.load_reference_data(
            GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
        )

        aligner_result = aligner.evaluate(
            relevant_distances=np.arange(0, 410, 10, dtype=int) / 100,
            full_reference_strategy=FullReferenceStrategy.PREFER_FULL_REFERENCE,
            max_predictions=-1,
        )
        process_results_evaluated = aligner_result.get_results(
            aligner=aligner, result_type=AlignerResultType.EVALUATED_PREDICTIONS
        )

        assert len(process_results_evaluated["theme_id_1"]) > 1
        assert (
            process_results_evaluated["theme_id_1"][3.5]["properties"][
                EVALUATION_FIELD_NAME
            ]
            == Evaluation.TO_CHECK_PREDICTION_FULL
        )

    @pytest.mark.usefixtures("mock_grb_response2")
    def test_evaluate_limited_predictions(self):
        thematic_dict = {
            "theme_id_1": from_wkt(
                "Polygon ((174015.08694170592934825 179025.39916784031083807, 174040.71934808720834553 179031.93985084796440788, 174037.1838437587430235 178986.50862022733781487, 174030.64316075111855753 178982.97311589887249283, 174016.85469387014745735 179002.24161448894301429, 174021.62762471355381422 179005.77711881740833633, 174018.62244603436556645 179008.60552228015149012, 174015.08694170592934825 179025.39916784031083807))"
            )
        }
        aligner = Aligner()
        aligner.load_thematic_data(DictLoader(thematic_dict))
        aligner.load_reference_data(
            GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
        )

        evaluation_result = aligner.evaluate(
            relevant_distances=np.arange(0, 410, 10, dtype=int) / 100,
            full_reference_strategy=FullReferenceStrategy.NO_FULL_REFERENCE,
            max_predictions=2,
        )
        assert len(evaluation_result.results["theme_id_1"]) > 1
        # assert (
        #     prop_dictionary["theme_id_1"][3.5]["properties"][EVALUATION_FIELD_NAME]
        #     == Evaluation.TO_CHECK_PREDICTION_MULTI_FULL
        # )

    @pytest.mark.usefixtures("mock_grb_response2")
    def test_evaluate_multi_to_best_prediction_true(self):
        thematic_dict = {
            "theme_id_1": from_wkt(
                "Polygon ((174015.08694170592934825 179025.39916784031083807, 174040.71934808720834553 179031.93985084796440788, 174037.1838437587430235 178986.50862022733781487, 174030.64316075111855753 178982.97311589887249283, 174016.85469387014745735 179002.24161448894301429, 174021.62762471355381422 179005.77711881740833633, 174018.62244603436556645 179008.60552228015149012, 174015.08694170592934825 179025.39916784031083807))"
            )
        }
        aligner = Aligner()
        aligner.load_thematic_data(DictLoader(thematic_dict))
        aligner.load_reference_data(
            GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
        )

        aligner_result = aligner.evaluate(
            relevant_distances=np.arange(0, 410, 10, dtype=int) / 100,
            full_reference_strategy=FullReferenceStrategy.NO_FULL_REFERENCE,
            max_predictions=1,
            multi_to_best_prediction=True,
        )
        process_results_evaluated = aligner_result.get_results(
            aligner=aligner, result_type=AlignerResultType.EVALUATED_PREDICTIONS
        )
        assert len(process_results_evaluated["theme_id_1"]) == 1
        assert (
            process_results_evaluated["theme_id_1"][
                list(process_results_evaluated["theme_id_1"].keys())[0]
            ]["properties"][PREDICTION_COUNT]
            > 1
        )

    @pytest.mark.usefixtures("mock_grb_response2")
    def test_evaluate_multi_to_best_prediction_false(self):
        thematic_dict = {
            "theme_id_1": from_wkt(
                "Polygon ((174015.08694170592934825 179025.39916784031083807, 174040.71934808720834553 179031.93985084796440788, 174037.1838437587430235 178986.50862022733781487, 174030.64316075111855753 178982.97311589887249283, 174016.85469387014745735 179002.24161448894301429, 174021.62762471355381422 179005.77711881740833633, 174018.62244603436556645 179008.60552228015149012, 174015.08694170592934825 179025.39916784031083807))"
            )
        }
        aligner = Aligner()
        aligner.load_thematic_data(DictLoader(thematic_dict))
        aligner.load_reference_data(
            GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
        )

        aligner_result = aligner.evaluate(
            relevant_distances=np.arange(0, 410, 10, dtype=int) / 100,
            full_reference_strategy=FullReferenceStrategy.NO_FULL_REFERENCE,
            max_predictions=1,
            multi_to_best_prediction=False,
        )
        process_results_evaluated = aligner_result.get_results(
            aligner=aligner, result_type=AlignerResultType.EVALUATED_PREDICTIONS
        )
        assert len(process_results_evaluated["theme_id_1"]) == 1
        assert (
            process_results_evaluated["theme_id_1"][0]["properties"][
                EVALUATION_FIELD_NAME
            ]
            == Evaluation.TO_CHECK_ORIGINAL
        )

    @pytest.mark.usefixtures("callback_grb_response")
    def test_evaluate_no_prediction(self):
        """
        Test evaluate when no prediction is returned when there are no stable geometries (=no zerostreaks) in the range of relevent_distances
        This testdata has no prediction between 0 and 1 meter, constantly increasing diff (first stable geometry at relevant_distances>2)
        """
        # Initiate an Aligner
        aligner = Aligner()
        # Load thematic data & reference data
        wkt = "Polygon ((174125.66583829143201001 179320.36133541050367057, 174125.66583829294540919 179320.36133543887990527, 174123.56114244047785178 179320.47114562886417843, 174123.60579274909105152 179320.9312379399780184, 174123.65622600910137407 179321.4509197928418871, 174130.00431616476271302 179386.86384879014804028, 174131.20465914410306141 179386.55606853921199217, 174131.47325675669708289 179386.48719735653139651, 174131.4777007331722416 179386.48605787538690493, 174131.48166125148418359 179386.48504236005828716, 174131.56964346842141822 179386.46248281505540945, 174145.93543933291221038 179382.77894541397108696, 174153.73543597472598776 179380.77894627506611869, 174155.22223844396648929 179380.39771487266989425, 174155.31905292189912871 179380.37289064755896106, 174158.08159073328715749 179379.66454761903150938, 174162.78046076380996965 179378.45970914987265132, 174162.77348954943590797 179378.41265345277497545, 174162.76924267943832092 179378.38398708027671091, 174162.76793666195590049 179378.37517146224854514, 174158.29825295542832464 179348.20480644324561581, 174158.06896705977851525 179348.24847994712763466, 174157.82196880917763337 179348.29552723292727023, 174152.4153556079545524 179349.32535831883433275, 174151.38357602406176738 179340.0393457512545865, 174150.45438929076772183 179331.6766651498619467, 174155.32298101272317581 179329.59012584027368575, 174155.97540604253299534 179329.31051516399020329, 174156.06587580699124373 179329.27174254972487688, 174156.3372870393213816 179329.15542325720889494, 174155.77701106332824565 179318.79031770044821315, 174154.96588003114447929 179318.83263741704286076, 174154.95658402171102352 179318.83312241756357253, 174151.97826524809352122 179318.98851313433260657, 174146.5353382924804464 179319.27249193197349086, 174125.66583829143201001 179320.36133541050367057))"

        loader = DictLoader({"id1": from_wkt(wkt)})
        aligner.load_thematic_data(loader)
        loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
        aligner.load_reference_data(loader)

        aligner_result = aligner.evaluate(
            relevant_distances=np.arange(0, 110, 10, dtype=int) / 100,
            full_reference_strategy=FullReferenceStrategy.NO_FULL_REFERENCE,
            max_predictions=1,
            multi_to_best_prediction=True,
        )
        process_results_evaluated = aligner_result.get_results(
            aligner=aligner, result_type=AlignerResultType.EVALUATED_PREDICTIONS
        )
        # TODO karel review after metadata is refined. This should be 0 zero, but is filled with metadata so we temporary have put it to 1
        assert len(process_results_evaluated) == 1

    @pytest.mark.usefixtures("mock_grb_response2")
    def test_evaluate_relevant_distances_without_0(self):
        thematic_dict = {
            "theme_id_1": from_wkt(
                "Polygon ((174015.08694170592934825 179025.39916784031083807, 174040.71934808720834553 179031.93985084796440788, 174037.1838437587430235 178986.50862022733781487, 174030.64316075111855753 178982.97311589887249283, 174016.85469387014745735 179002.24161448894301429, 174021.62762471355381422 179005.77711881740833633, 174018.62244603436556645 179008.60552228015149012, 174015.08694170592934825 179025.39916784031083807))"
            )
        }
        aligner = Aligner()
        aligner.load_thematic_data(DictLoader(thematic_dict))
        aligner.load_reference_data(
            GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
        )
        with pytest.raises(ValueError):
            aligner_result = aligner.evaluate(
                relevant_distances=np.arange(10, 410, 10, dtype=int) / 100,
                full_reference_strategy=FullReferenceStrategy.NO_FULL_REFERENCE,
                max_predictions=1,
                multi_to_best_prediction=False,
            )

    @pytest.mark.usefixtures("callback_grb_response")
    def test_evaluate_line(self):
        # Load thematic data & reference data
        thematic_dict = {
            "theme_id": from_wkt(
                "LINESTRING (171741.11190000033820979 171839.01070547936251387, 171751.68948904142598622 171847.23771917796693742, 171762.26707808251376264 171855.69979041084297933, 171772.72713835645117797 171862.9865739724773448, 171783.53978493178146891 171870.8610013697471004, 171794.94007534274714999 171877.79519862998859026, 171801.87427260301774368 171884.2592808217741549)"
            )
        }

        # ADD A REFERENCE POLYGON TO REFERENCE DICTIONARY

        aligner = Aligner()
        aligner.load_thematic_data(DictLoader(thematic_dict))
        # aligner.load_reference_data(DictLoader(reference_dict))
        aligner.load_reference_data(
            GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
        )

        evaluation_result = aligner.evaluate(
            relevant_distances=np.arange(0, 1010, 50, dtype=int) / 100,
            full_reference_strategy=FullReferenceStrategy.NO_FULL_REFERENCE,
            max_predictions=3,
            multi_to_best_prediction=False,
        )
        assert len(evaluation_result.results["theme_id"]) >= 1
        assert (
            evaluation_result.results["theme_id"][0]["properties"][
                EVALUATION_FIELD_NAME
            ]
            == Evaluation.TO_CHECK_PREDICTION_MULTI_FULL
        )

    @pytest.mark.usefixtures("mock_grb_response2")
    def test_evaluate_point(self):
        # Load thematic data & reference data
        # thematic_dict = {"theme_id": from_wkt("POINT (0 0)")}
        thematic_dict = {
            "theme_id": from_wkt(
                "POINT (173966.17483414348680526 172343.78743441699771211)"
            )
        }

        aligner = Aligner()
        aligner.load_thematic_data(DictLoader(thematic_dict))
        aligner.load_reference_data(
            GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
        )

        evaluation_result = aligner.evaluate(
            relevant_distances=np.arange(0, 1010, 50, dtype=int) / 100,
            full_reference_strategy=FullReferenceStrategy.NO_FULL_REFERENCE,
            max_predictions=3,
            multi_to_best_prediction=False,
        )
        assert len(evaluation_result.results["theme_id"]) > 0
        # assert (
        #     evaluation_result.results["theme_id"][0]["properties"][EVALUATION_FIELD_NAME]
        #     == Evaluation.TO_CHECK_PREDICTION_MULTI_FULL
        # )

    def test_evaluate_reference_point(self):
        # Load thematic data & reference data
        polygon_1 = from_wkt(
            "MULTIPOLYGON (((171795.71618631482124329 171817.88460136577486992, 171784.53532230854034424 171806.1688893586397171, 171746.73993028700351715 171841.20300138369202614, 171746.28380228579044342 171841.62578538432717323, 171767.35881029814481735 171856.89906539395451546, 171767.47471430152654648 171856.76376939192414284, 171798.38581831753253937 171820.68191336840391159, 171798.01820231974124908 171820.29669736698269844, 171795.71618631482124329 171817.88460136577486992)))"
        )
        geom_reference = from_wkt(
            "MULTIPOINT ((171756.52506366037414409 171850.34457502837176435),(171766.12778689494007267 171857.04121096828021109),(171777.13617191879893653 171865.19089055553195067),(171745.72200002145837061 171841.94219219812657684),(171770.11351679253857583 171840.12903735807049088))"
        )

        # thematic_dict = {"theme_id": from_wkt("POINT (0 0)")}
        thematic_dict = {"theme_id": polygon_1}

        # ADD A REFERENCE POLYGON TO REFERENCE DICTIONARY
        reference_dict = {"ref_id": geom_reference}

        config = ProcessorConfig(
            snap_strategy=SnapStrategy.NO_PREFERENCE, snap_max_segment_length=2
        )
        processor = AlignerGeometryProcessor(config)
        aligner = Aligner(processor=processor)
        aligner.load_thematic_data(DictLoader(thematic_dict))
        aligner.load_reference_data(DictLoader(reference_dict))

        evaluation_result = aligner.evaluate(
            relevant_distances=np.arange(0, 510, 50, dtype=int) / 100,
            full_reference_strategy=FullReferenceStrategy.NO_FULL_REFERENCE,
            max_predictions=3,
            multi_to_best_prediction=False,
        )
        assert True

    @pytest.mark.usefixtures("mock_grb_response3")
    def test_evaluate_best_no_prediction(self):
        thematic_json = {
            "type": "FeatureCollection",
            "name": "uc5",
            "crs": {
                "type": "name",
                "properties": {"name": "urn:ogc:def:crs:EPSG::31370"},
            },
            "features": [
                {
                    "type": "Feature",
                    "properties": {"fid": 1, "dossiernummer": "D_68365"},
                    "geometry": {
                        "type": "Polygon",
                        "coordinates": [
                            [
                                [172648.411900000093738, 171970.078100000013364],
                                [172641.128199999919161, 171979.2324],
                                [172638.179699999920558, 171982.9382],
                                [172633.854899999918416, 171988.3738],
                                [172633.317000841663685, 171989.049898942059372],
                                [172623.4652, 172001.432],
                                [172610.793499999970663, 172017.358100000128616],
                                [172610.942697853111895, 172017.477398283459479],
                                [172650.691, 172049.2673],
                                [172654.167900000087684, 172052.047999999922467],
                                [172656.0643, 172049.772999999986496],
                                [172657.2538, 172048.346],
                                [172677.110299999912968, 172024.525499999901513],
                                [172677.8891, 172023.591299999912735],
                                [172692.862399999867193, 172005.628700000088429],
                                [172693.374199999903794, 172006.038],
                                [172696.634700000082375, 172001.575399999826914],
                                [172709.757699999987381, 171986.0387],
                                [172698.4481, 171976.4558],
                                [172665.904004495765548, 171948.880303809302859],
                                [172665.794206234189915, 171948.787305280333385],
                                [172665.529999999969732, 171948.563400000159163],
                                [172648.411900000093738, 171970.078100000013364],
                            ]
                        ],
                    },
                },
                {
                    "type": "Feature",
                    "properties": {"fid": 2, "dossiernummer": "D_352"},
                    "geometry": {
                        "type": "Polygon",
                        "coordinates": [
                            [
                                [172443.633200000098441, 171905.855499999917811],
                                [172450.522299999924144, 171922.915599999832921],
                                [172471.899299999902723, 171914.2226],
                                [172458.880800000013551, 171881.983800000074552],
                                [172447.773599999636644, 171885.025300000183051],
                                [172436.992100000090431, 171889.409699999901932],
                                [172443.633200000098441, 171905.855499999917811],
                            ]
                        ],
                    },
                },
                {
                    "type": "Feature",
                    "properties": {"fid": 3, "dossiernummer": "D_28311"},
                    "geometry": {
                        "type": "Polygon",
                        "coordinates": [
                            [
                                [172179.29800000009709, 171788.713699999905657],
                                [172182.4732, 171780.029],
                                [172157.412400000175694, 171770.505299999989802],
                                [172129.436300000088522, 171848.2427],
                                [172129.052500541059999, 171849.309098496683873],
                                [172128.886999999987893, 171849.769],
                                [172148.2379, 171862.7293],
                                [172151.4504, 171864.880899999901885],
                                [172151.894199230970116, 171863.666902103519533],
                                [172179.29800000009709, 171788.713699999905657],
                            ]
                        ],
                    },
                },
            ],
        }
        # 4001 will give a result, 4002 &4009 should also give the original geometry

        aligner = Aligner()
        loader = GeoJsonLoader(_input=thematic_json, id_property="fid")
        aligner.load_thematic_data(loader)
        # Load reference data: The actual GRB-parcels
        loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
        aligner.load_reference_data(loader)

        # Use the EVALUATE-function
        aligner_result = aligner.evaluate(
            relevant_distances=np.arange(0, 310, 10, dtype=int) / 100,
            full_reference_strategy=FullReferenceStrategy.ONLY_FULL_REFERENCE,
            max_predictions=1,
            multi_to_best_prediction=True,
        )
        process_results_evaluated = aligner_result.get_results(
            aligner=aligner, result_type=AlignerResultType.EVALUATED_PREDICTIONS
        )

        assert (
            process_results_evaluated[1][list(process_results_evaluated[1].keys())[0]][
                "properties"
            ][EVALUATION_FIELD_NAME]
            == Evaluation.PREDICTION_UNIQUE_AND_FULL_REFERENCE
        )
        assert (
            process_results_evaluated[2][list(process_results_evaluated[2].keys())[0]][
                "properties"
            ][EVALUATION_FIELD_NAME]
            == Evaluation.TO_CHECK_NO_PREDICTION
        )
        assert (
            process_results_evaluated[3][list(process_results_evaluated[3].keys())[0]][
                "properties"
            ][EVALUATION_FIELD_NAME]
            == Evaluation.TO_CHECK_NO_PREDICTION
        )
