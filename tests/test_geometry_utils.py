import unittest
from shapely.geometry import Point, Polygon
from brdr.geometry_utils import (
    buffer_neg_pos,
    buffer_neg,
    buffer_pos,
    safe_union,
    safe_intersection,
    safe_difference,
    safe_symmetric_difference,
    grid_bounds,
)


class TestBuffering(unittest.TestCase):
    def test_buffer_neg_pos_point(self):
        """Tests buffer_neg_pos with a point geometry."""
        point = Point(0, 0)
        result = buffer_neg_pos(point, 1.0)
        self.assertEqual(result, Polygon())

    def test_buffer_neg_pos_polygon(self):
        """Tests buffer_neg_pos with a polygon geometry."""
        polygon = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
        result = buffer_neg_pos(polygon, 0.5)
        self.assertTrue(result.is_valid)

    def test_buffer_neg(self):
        """Tests buffer_neg with a point geometry."""
        point = Point(0, 0)
        result = buffer_neg(point, 1.0)
        self.assertEqual(result, Polygon())

    def test_buffer_pos(self):
        """Tests buffer_pos with a point geometry."""
        polygon = Polygon(((-1, -1), (-1, 1), (1, 1), (1, -1)))
        result = buffer_pos(polygon, 1.0)
        self.assertEqual(result, Polygon(((-2, -2), (-2, 2), (2, 2), (2, -2), (-2, -2))))


class TestSafeOperations(unittest.TestCase):
    #TODO: add cases where GEOS-exception occurs
    def test_safe_union(self):
        """Tests safe_union with two points."""
        polygon_a = Polygon(((0, 0), (1, 0), (1, 1), (0, 1), (0, 0)))
        polygon_b = Polygon(((1, 0), (2, 0), (2, 1), (1, 1), (1, 0)))
        result = safe_union(polygon_a, polygon_b)
        print(result)
        self.assertEqual(result.area, 2)

    def test_safe_union_geos_exception(self):
        """Tests safe_union with two points where GEOS exception occurs."""
        #TODO search for this case
        polygon_a = Polygon(((0, 0), (1, 0), (1, 1), (0, 1), (0, 0)))
        polygon_b = Polygon(((1, 0), (2, 0), (2, 1), (1, 1), (1, 0)))
        result = safe_union(polygon_a, polygon_b)
        print(result)
        self.assertEqual(result.area, 2)

    def test_safe_intersection(self):
        """Tests safe_intersection with two disjoint polygons."""
        polygon_a = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
        polygon_b = Polygon([(2, 2), (3, 2), (3, 3), (2, 3)])
        result = safe_intersection(polygon_a, polygon_b)
        self.assertEqual(result, Polygon())

    def test_safe_intersection_geos_exception(self):
        """Tests safe_intersection with two points where GEOS exception occurs."""
        #TODO search for this case
        polygon_a = Polygon(((0, 0), (1, 0), (1, 1), (0, 1), (0, 0)))
        polygon_b = Polygon(((1, 0), (2, 0), (2, 1), (1, 1), (1, 0)))
        result = safe_intersection(polygon_a, polygon_b)
        print(result)
        self.assertEqual(result.area, 0)
    def test_safe_difference(self):
        """Tests safe_difference with a polygon contained within another."""
        polygon_a = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
        polygon_b = Polygon([(0.2, 0.2), (0.8, 0.2), (0.8, 0.8), (0.2, 0.8)])
        result = safe_difference(polygon_a, polygon_b)
        self.assertTrue(result.is_valid)

    def test_safe_difference_geos_exception(self):
        """Tests safe_difference with a polygon contained within another."""
        # TODO search for this case
        polygon_a = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
        polygon_b = Polygon([(0.2, 0.2), (0.8, 0.2), (0.8, 0.8), (0.2, 0.8)])
        result = safe_difference(polygon_a, polygon_b)
        self.assertTrue(result.is_valid)

    def test_safe_symmetric_difference(self):
        """Tests safe_symmetric_difference with two disjoint polygons."""
        polygon_a = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
        polygon_b = Polygon([(2, 2), (3, 2), (3, 3), (2, 3)])
        result = safe_symmetric_difference(polygon_a, polygon_b)
        self.assertTrue(result.is_valid)

    def test_safe_symmetric_difference_geos_exception(self):
        # TODO search for this case
        """Tests safe_symmetric_difference with two disjoint polygons."""
        polygon_a = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
        polygon_b = Polygon([(2, 2), (3, 2), (3, 3), (2, 3)])
        result = safe_symmetric_difference(polygon_a, polygon_b)
        self.assertTrue(result.is_valid)


class TestGridBounds(unittest.TestCase):
    def test_grid_bounds_empty_polygon(self):
        """Tests grid_bounds with an empty polygon."""
        polygon = Polygon()
        result = grid_bounds(polygon, 1.0)
        self.assertEqual(result, polygon)

    def test_grid_bounds_small_grid(self):
        """Tests grid_bounds with a small area not requiring grid division."""
        polygon = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
        result = grid_bounds(polygon, 2.0)
        self.assertEqual(len(result), 1)

    def test_grid_bounds_grid_division(self):
        """Tests grid_bounds with an area requiring grid"""
        polygon = Polygon([(0, 0), (0, 5), (5, 5), (5, 0), (0, 0)])
        result = grid_bounds(polygon, 1.0)
        self.assertEqual(len(result), 25)


