import json
import unittest

import numpy as np
import pytest
from shapely import to_geojson
from shapely.geometry import shape

from brdr.aligner import Aligner
from brdr.be.grb.enums import GRBType
from brdr.be.grb.loader import GRBActualLoader
from brdr.configs import ProcessorConfig
from brdr.enums import DiffMetric
from brdr.enums import OpenDomainStrategy
from brdr.loader import DictLoader
from brdr.processor import AlignerGeometryProcessor


class TestExamples(unittest.TestCase):

    @pytest.mark.usefixtures("callback_grb_response")
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

        config = ProcessorConfig(
            area_limit=100000, od_strategy=OpenDomainStrategy.SNAP_INNER_SIDE
        )
        processor = AlignerGeometryProcessor(config)
        aligner = Aligner(processor=processor)

        geometry = shape(contour)

        aligner.load_thematic_data(DictLoader({"input_id": geometry}))
        aligner.load_reference_data(
            GRBActualLoader(
                grb_type=referentielaag_type, partition=1000, aligner=aligner
            )
        )

        series = np.arange(0, 61, 1, dtype=float) / 10

        process_result = aligner.process(relevant_distances=series)
        dict_diffs = aligner.get_difference_metrics_for_thematic_data(
            process_result.results,
            aligner.thematic_data,
            DiffMetric.SYMMETRICAL_AREA_CHANGE,
        )

        dict_diffs = dict_diffs["input_id"]
        serial_dict = {}
        dict_results = process_result.results["input_id"]
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
