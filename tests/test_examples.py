import unittest

import numpy as np

from brdr.aligner import Aligner
from brdr.enums import GRBType
from brdr.grb import GRBActualLoader
from brdr.loader import DictLoader
from brdr.loader import GeoJsonLoader
from brdr.oe import get_oe_dict_by_ids, OnroerendErfgoedLoader
from brdr.utils import diffs_from_dict_series
from brdr.utils import get_breakpoints_zerostreak


class TestExamples(unittest.TestCase):

    def test_example_131635(self):
        # EXAMPLE
        aligner = Aligner()
        dict_theme = get_oe_dict_by_ids([131635])
        loader = DictLoader(dict_theme)
        aligner.load_thematic_data(loader)
        loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
        aligner.load_reference_data(loader)
        rel_dist = 2
        aligner.process(relevant_distance=rel_dist, od_strategy=4)

    def test_example_combined_borders_adp_gbg(self):
        aligner = Aligner()
        loader = OnroerendErfgoedLoader([131635])
        aligner.load_thematic_data(loader)
        adp_loader = GRBActualLoader(
            grb_type=GRBType.ADP, partition=1000, aligner=aligner
        )
        gbg_loader = GRBActualLoader(
            grb_type=GRBType.GBG, partition=1000, aligner=aligner
        )
        dict_ref, dict_ref_properties_adp, source_adp = adp_loader.load_data()
        dict_ref2, dict_ref_properties_gbg, source_gbg = gbg_loader.load_data()
        dict_ref.update(dict_ref2)  # combine 2 dictionaries
        # make a polygonized version of the reference data with non-overlapping polygons
        aligner.load_reference_data(DictLoader(dict_ref))

        rel_dist = 2
        result_dict = aligner.process(relevant_distance=rel_dist, od_strategy=4)
        for process_results in result_dict.values():
            aligner.get_brdr_formula(process_results[rel_dist]["result"])

    def test_example_multipolygon(self):
        aligner0 = Aligner()
        testdata = {
            "type": "FeatureCollection",
            "name": "themelayer",
            "crs": {
                "type": "name",
                "properties": {"name": "urn:ogc:def:crs:EPSG::31370"},
            },
            "features": [
                {
                    "type": "Feature",
                    "properties": {"id": 600, "theme_identifier": "600"},
                    "geometry": {
                        "type": "MultiPolygon",
                        "coordinates": [
                            [
                                [
                                    [172907.478557254769839, 171306.329770992975682],
                                    [172907.423273425694788, 171307.187050310545601],
                                    [172907.258293797436636, 171308.030809821182629],
                                    [172906.986220197344664, 171308.847742933343397],
                                    [172906.611343388562091, 171309.62496612383984],
                                    [172906.139575402339688, 171310.350222118664533],
                                    [172905.578356301615713, 171311.012073197518475],
                                    [172904.936536846886156, 171311.600081573560601],
                                    [172904.224238914379384, 171312.104974003392272],
                                    [172903.452695867919829, 171312.518788031826261],
                                    [172902.634075402340386, 171312.834997564437799],
                                    [172901.781287651334424, 171313.048615787964081],
                                    [172900.907781587186037, 171313.156273815431632],
                                    [172900.027332922356436, 171313.156273815431632],
                                    [172899.153826858208049, 171313.048615787964081],
                                    [172898.301039107202087, 171312.834997564437799],
                                    [172897.482418641622644, 171312.518788031826261],
                                    [172896.710875595163088, 171312.104974003392272],
                                    [172895.998577662656317, 171311.600081573560601],
                                    [172895.356758207926759, 171311.012073197518475],
                                    [172894.795539107202785, 171310.350222118664533],
                                    [172894.323771120980382, 171309.62496612383984],
                                    [172893.948894312197808, 171308.847742933343397],
                                    [172893.676820712105837, 171308.030809821182629],
                                    [172893.511841083847685, 171307.187050310545601],
                                    [172893.456557254772633, 171306.329770992975682],
                                    [172893.511841083847685, 171305.472491675405763],
                                    [172893.676820712105837, 171304.628732164768735],
                                    [172893.948894312197808, 171303.811799052607967],
                                    [172894.323771120980382, 171303.034575862111524],
                                    [172894.795539107202785, 171302.309319867286831],
                                    [172895.356758207926759, 171301.647468788432889],
                                    [172895.998577662656317, 171301.059460412390763],
                                    [172896.710875595163088, 171300.554567982559092],
                                    [172897.482418641622644, 171300.140753954125103],
                                    [172898.301039107202087, 171299.824544421513565],
                                    [172899.153826858208049, 171299.610926197987283],
                                    [172900.027332922356436, 171299.503268170519732],
                                    [172900.907781587186037, 171299.503268170519732],
                                    [172901.781287651334424, 171299.610926197987283],
                                    [172902.634075402340386, 171299.824544421513565],
                                    [172903.452695867919829, 171300.140753954125103],
                                    [172904.224238914379384, 171300.554567982559092],
                                    [172904.936536846886156, 171301.059460412390763],
                                    [172905.578356301615713, 171301.647468788432889],
                                    [172906.139575402339688, 171302.309319867286831],
                                    [172906.611343388562091, 171303.034575862111524],
                                    [172906.986220197344664, 171303.811799052607967],
                                    [172907.258293797436636, 171304.628732164768735],
                                    [172907.423273425694788, 171305.472491675405763],
                                    [172907.478557254769839, 171306.329770992975682],
                                ]
                            ],
                            [
                                [
                                    [172859.258396949415328, 171379.685496183810756],
                                    [172893.114885499031516, 171367.032061069301562],
                                    [172879.777480918884976, 171334.201526718155947],
                                    [172847.630916033376707, 171347.025954199081752],
                                    [172859.258396949415328, 171379.685496183810756],
                                ]
                            ],
                            [
                                [
                                    [172908.846183208952425, 171363.954198473889846],
                                    [172929.536259544838686, 171356.259541985346004],
                                    [172926.116412216593744, 171336.766412214346929],
                                    [172902.177480918879155, 171346.341984733415302],
                                    [172908.846183208952425, 171363.954198473889846],
                                ]
                            ],
                            [
                                [
                                    [172911.411068705143407, 171403.96641221435857],
                                    [172908.846183208952425, 171402.598473283054773],
                                    [172908.333206109731691, 171398.323664122755872],
                                    [172910.214122140256222, 171395.75877862656489],
                                    [172915.172900766221574, 171395.929770992981503],
                                    [172918.592748094466515, 171399.007633588393219],
                                    [172917.908778628800064, 171401.914503817417426],
                                    [172915.685877865442308, 171404.308396947191795],
                                    [172912.608015270030592, 171404.992366412829142],
                                    [172911.411068705143407, 171403.96641221435857],
                                ]
                            ],
                        ],
                    },
                }
            ],
        }

        # Load thematic data
        aligner0.load_thematic_data(
            GeoJsonLoader(_input=testdata, id_property="theme_identifier")
        )
        aligner0.multi_as_single_modus = True
        dict_thematic = aligner0.dict_thematic
        aligner0.load_thematic_data(DictLoader(dict_thematic))

        # gebruik de actuele adp-percelen adp= administratieve percelen
        aligner = Aligner()
        aligner.load_thematic_data(DictLoader(dict_thematic))
        aligner.load_reference_data(
            GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
        )

        _, dict_predictions, _ = aligner.predictor()

        self.assertGreater(len(dict_predictions), 0)
        fcs = aligner.get_results_as_geojson(formula=True)
        self.assertEqual(len(fcs), 6)

    def test_example_wanted_changes(self):
        aligner = Aligner()
        # Load thematic data & reference data
        dict_theme = get_oe_dict_by_ids([131635])
        aligner.load_thematic_data(DictLoader(dict_theme))
        aligner.load_reference_data(
            GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
        )

        # Example how to use the Aligner
        rel_dist = 2
        aligner.process(relevant_distance=rel_dist, od_strategy=4)

        # Example how to use a series (for histogram)
        series = np.arange(0, 300, 10, dtype=int) / 100
        dict_series = aligner.process(series, 4, 50)
        resulting_areas = diffs_from_dict_series(dict_series, aligner.dict_thematic)
        for key in resulting_areas:
            if len(resulting_areas[key]) == len(series):
                lst_diffs = list(resulting_areas[key].values())
                extremes, zero_streak = get_breakpoints_zerostreak(series, lst_diffs)
                print(str(key))
                for extremum in extremes:
                    print(f"{extremum[0]:.2f}, {extremum[1]:.2f} ({extremum[2]})")
                for st in zero_streak:
                    print(
                        f"{st[0]:.2f} - {st[1]:.2f} -{st[2]:.2f} - {st[3]:.2f}"
                        f" - startextreme {st[4]:.2f} "
                    )
                    aligner.process(relevant_distance=st[0], od_strategy=4)

    def test_example_predictor(self):
        aligner = Aligner()
        # Load thematic data & reference data
        dict_theme = get_oe_dict_by_ids([131635])
        aligner.load_thematic_data(DictLoader(dict_theme))
        aligner.load_reference_data(
            GRBActualLoader(aligner=aligner, grb_type=GRBType.GBG, partition=1000)
        )
        series = np.arange(0, 300, 10, dtype=int) / 100
        # predict which relevant distances are interesting to propose as resulting
        # geometry

        _, dict_predictions, _ = aligner.predictor(
            relevant_distances=series, od_strategy=4, threshold_overlap_percentage=50
        )
        for key in dict_predictions.keys():
            assert key in dict_predictions.keys()
            continue
