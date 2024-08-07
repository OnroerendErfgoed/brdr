import unittest

from shapely import Polygon

from brdr.aligner import Aligner
from brdr.enums import GRBType
from brdr.loader import DictLoader
from brdr.loader import GRBActualLoader
from brdr.utils import (
    get_oe_dict_by_ids,
)


class TestExamples(unittest.TestCase):

    def test_load_data(self):
        # EXAMPLE
        aligner = Aligner()

        dict_theme = get_oe_dict_by_ids([131635])
        thematic_loader = DictLoader(data_dict=dict_theme)
        reference_loader = GRBActualLoader(
            grb_type=GRBType.ADP, aligner=aligner, partition=0
        )

        aligner.dict_thematic = thematic_loader.load_data()
        old_thematic_data = aligner.dict_thematic
        aligner.dict_reference = reference_loader.load_data()
        assert aligner.dict_reference is not None

    def test_partition(self):
        # Test partition function
        delta = 2.0
        sample_geom = Polygon([(0, 0), (0, 10), (10, 10), (10, 0)])
        filtered_partitions = GRBActualLoader.partition(sample_geom, delta)

        # Check if the result is a list of Polygon objects
        self.assertIsInstance(filtered_partitions, list)
        for partition in filtered_partitions:
            self.assertIsInstance(partition, Polygon)
