import unittest

import numpy as np
from shapely import from_wkt

from brdr.aligner import Aligner
from brdr.configs import ProcessorConfig
from brdr.loader import DictLoader
from brdr.processor import NetworkGeometryProcessor


class TestProcessor(unittest.TestCase):
    def setUp(self):
        pass

    def test_networkgeometryprocessor(self):
        # Load thematic data & reference data
        thematic_dict = {"theme_id": from_wkt("POLYGON ((0 0, 0 9, 5 10, 10 0, 0 0))")}
        # ADD A REFERENCE POLYGON TO REFERENCE DICTIONARY
        reference_dict = {"ref_id": from_wkt("POLYGON ((0 1, 0 10,8 10,10 1,0 1))")}
        # LOAD THEMATIC DICTIONARY
        processor = NetworkGeometryProcessor(config=ProcessorConfig())
        aligner = Aligner(processor=processor)
        aligner.load_thematic_data(DictLoader(thematic_dict))
        # LOAD REFERENCE DICTIONARY
        aligner.load_reference_data(DictLoader(reference_dict))
        series = np.arange(0, 310, 10, dtype=int) / 100
        # predict which relevant distances are interesting to propose as resulting
        # geometry

        prediction_result = aligner.predict(series)
        assert len(prediction_result.results) == len(thematic_dict)
