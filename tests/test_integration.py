import unittest

import numpy as np

from brdr.aligner import Aligner
from brdr.enums import GRBType
from brdr.grb import GRBActualLoader
from brdr.loader import DictLoader, GeoJsonLoader
from brdr.utils import diffs_from_dict_series
from brdr.utils import get_breakpoints_zerostreak
from brdr.utils import get_oe_dict_by_ids
from brdr.utils import multipolygons_to_singles

import json

import numpy as np
from brdr.aligner import Aligner
from brdr.enums import DiffMetric
from brdr.enums import GRBType
from brdr.enums import OpenbaarDomeinStrategy
from brdr.grb import GRBActualLoader
from brdr.loader import DictLoader
from brdr.typings import GeoJSONGeometry
from brdr.utils import diffs_from_dict_series
from shapely import to_geojson
from shapely.geometry import shape


class TestExamples(unittest.TestCase):



    def test_webservice_brdr(self
    ):
        """
        Code used in brdr-webservice
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
        openbaardomein_strategy = OpenbaarDomeinStrategy.SNAP_SINGLE_SIDE

        aligner = Aligner(
            area_limit=100000
        )

        geometry = shape(contour)

        aligner.load_thematic_data(DictLoader({"input": geometry}))
        aligner.load_reference_data(
            GRBActualLoader(grb_type=referentielaag_type, partition=1000, aligner=aligner)
        )

        series = np.arange(0, 61, 1, dtype=float) / 10

        dict_series = aligner.process_series(series, openbaardomein_strategy, 50)
        dict_diffs = diffs_from_dict_series(
            dict_series, aligner.dict_thematic, DiffMetric.CHANGES_AREA
        )

        dict_diffs = dict_diffs["input"]
        dict_series = {
            rel_dist: {
                "result": json.loads(to_geojson(dict_results["input"]["result"])),
                "result_diff_min": json.loads(
                    to_geojson(dict_results["input"]["result_diff_min"])
                ),
                "result_diff_plus": json.loads(
                    to_geojson(dict_results["input"]["result_diff_plus"])
                ),
            }
            for rel_dist, dict_results in dict_series.items()
        }

        return {
            "series": dict_series,
            "diffs": dict_diffs,
        }