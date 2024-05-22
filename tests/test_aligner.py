import unittest

from shapely import from_wkt
from shapely.geometry import Polygon

from brdr.aligner import Aligner
from brdr.geometry_utils import buffer_neg_pos
from brdr.geometry_utils import grid_bounds


class TestAligner(unittest.TestCase):
    def setUp(self):
        # Create a sample geometry for testing
        self.sample_aligner = Aligner()
        self.sample_geom = Polygon([(0, 0), (0, 10), (10, 10), (10, 0)])

    def test_buffer_neg_pos(self):
        # Check if the result of the _buffer_neg_pos gives an equal geometry
        out = buffer_neg_pos(self.sample_geom, 3)
        self.assertEqual(self.sample_geom, out)

    def test_grid_bounds_1(self):
        # Test _grid_bounds function
        delta = 1.0
        grid_partitions = grid_bounds(self.sample_geom, delta)

        # Check if the result is a list of Polygon objects
        self.assertIsInstance(grid_partitions, list)
        for partition in grid_partitions:
            self.assertIsInstance(partition, Polygon)

        # Add more specific tests based on your requirements

    def test_grid_bounds_2(self):
        # Test _grid_bounds function
        delta = 2.0
        grid_partitions = grid_bounds(self.sample_geom, delta)

        # Check if the result is a list of Polygon objects
        self.assertIsInstance(grid_partitions, list)
        for partition in grid_partitions:
            self.assertIsInstance(partition, Polygon)

    def test_get_last_version_date(self):
        # Check if the result of the _buffer_neg_pos gives an equal geometry
        Polygon([(0, 0), (0, 10), (10, 10), (10, 0)])
        geom = Polygon([(0, 0), (0, 10), (10, 10), (10, 0)])
        out = self.sample_aligner.get_last_version_date(geom)
        self.assertIsNone(out)
        geom = Polygon(
            [(170000, 170000), (170000, 170100), (170100, 170100), (170100, 170000)]
        )
        out = self.sample_aligner.get_last_version_date(geom)
        self.assertRegex(
            out, "^(?:19|20)\\d\\d-(?:0[1-9]|1[0-2])-(?:0[1-9]|[12][0-9]|3[01])$"
        )

    def test_partition(self):
        # Test partition function
        delta = 2.0
        filtered_partitions = self.sample_aligner.partition(
            self.sample_geom, delta
        )

        # Check if the result is a list of Polygon objects
        self.assertIsInstance(filtered_partitions, list)
        for partition in filtered_partitions:
            self.assertIsInstance(partition, Polygon)

        # Add more specific tests based on your requirements

    def test_get_formula_full_intersection(self):
        # Test when intersection equals reference geometry
        key = "a"
        ref_dict = {key: self.sample_geom}
        self.sample_aligner.load_reference_data_dict(ref_dict)
        res = self.sample_aligner.get_formula(self.sample_geom, with_geom=True)
        result = res[key]
        self.assertTrue(result["full"])
        self.assertEqual(result["percentage"], 100)

    def test_get_formula_partial_intersection(self):
        # Test when intersection is partial
        key = "a"
        ref_dict = {key: self.sample_geom.buffer(0.5)}
        self.sample_aligner.load_reference_data_dict(ref_dict)
        res = self.sample_aligner.get_formula(self.sample_geom, with_geom=True)
        result = res[key]
        self.assertFalse(result["full"])
        self.assertGreater(result["percentage"], 0)
        self.assertLess(result["percentage"], 100)

    def test_process_geometry(self):
        # Test if processed geometry is equal to reference geometry
        key_ref = "a"
        ref_dict = {key_ref: self.sample_geom}
        self.sample_aligner.load_reference_data_dict(ref_dict)
        r, rd, rdp, rdm, si, sd = self.sample_aligner.process_geometry(
            self.sample_geom.buffer(0.5)
        )
        self.assertTrue(from_wkt(r.wkt).equals(from_wkt(self.sample_geom.wkt)))
        self.assertFalse(rd.is_empty)
