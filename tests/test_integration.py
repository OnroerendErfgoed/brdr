import json
import unittest

import numpy as np
from shapely import to_geojson
from shapely.geometry import shape

from brdr.aligner import Aligner
from brdr.enums import DiffMetric
from brdr.enums import GRBType
from brdr.enums import OpenDomainStrategy
from brdr.grb import GRBActualLoader
from brdr.loader import DictLoader


class TestExamples(unittest.TestCase):

    def test_webservice_inventaris_brdr_integration(self):
        """
        Code used in inventaris-webservice for brdr-integration
        """
        contour = {
            "type": "MultiPolygon",
            "coordinates": [
                [
                    [
                        [174647.924166895216331, 170964.980246312363306],
                        [174693.204879119351972, 170943.531487890402786],
                        [174696.382472959638108, 170930.82111252922914],
                        [174678.905706838035258, 170901.428369506524177],
                        [174660.634542256389977, 170861.708446502889274],
                        [174641.56897921464406, 170820.399726579111302],
                        [174593.905071610264713, 170839.46528962085722],
                        [174614.559431572153699, 170881.568408004706725],
                        [174628.064205393398879, 170926.054721768799936],
                        [174647.924166895216331, 170964.980246312363306],
                    ]
                ]
            ],
        }
        referentielaag_type = GRBType.ADP
        openbaardomein_strategy = OpenDomainStrategy.SNAP_INNER_SIDE

        aligner = Aligner(area_limit=100000)

        geometry = shape(contour)

        aligner.load_thematic_data(DictLoader({"input_id": geometry}))
        aligner.load_reference_data(
            GRBActualLoader(
                grb_type=referentielaag_type, partition=1000, aligner=aligner
            )
        )

        series = np.arange(0, 61, 1, dtype=float) / 10

        dict_series = aligner.process(
            relevant_distances=series,
            od_strategy=openbaardomein_strategy,
            threshold_overlap_percentage=50,
        )
        dict_diffs = aligner.get_diff_metrics(
            dict_series, aligner.dict_thematic, DiffMetric.CHANGES_AREA
        )

        dict_diffs = dict_diffs["input_id"]
        serial_dict = {}
        dict_results = dict_series["input_id"]
        for rel_dist, process_results in dict_results.items():
            serial_dict[rel_dist] = {
                "result": json.loads(to_geojson(dict_results[rel_dist]["result"])),
                "result_diff_min": json.loads(
                    to_geojson(dict_results[rel_dist]["result_diff_min"])
                ),
                "result_diff_plus": json.loads(
                    to_geojson(dict_results[rel_dist]["result_diff_plus"])
                ),
            }
        assert len(serial_dict) == len(series)
        assert len(dict_diffs) == len(series)
