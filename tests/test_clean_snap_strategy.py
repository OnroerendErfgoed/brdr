import unittest

import networkx as nx
from shapely.geometry import LineString, MultiLineString

from brdr.enums import SnapStrategy
from brdr.geometry_utils import (
    bridge_with_straight_line,
    clean_pseudo_nodes_by_snap_strategy,
    find_best_path_in_network,
    get_end_coords,
)


class TestCleanSnapStrategy(unittest.TestCase):
    def test_get_end_coords_returns_original_reference_ends(self):
        original_reference = LineString([(0.0, 0.0), (2.0, 0.0), (4.0, 0.0)])
        ref_split = MultiLineString([[(0.0, 0.0), (2.0, 0.0)], [(2.0, 0.0), (4.0, 0.0)]])

        end_coords = get_end_coords(original_reference, ref_split)

        self.assertEqual(end_coords, {(0.0, 0.0), (4.0, 0.0)})

    def test_bridge_preserves_polyline_shape_when_removing_node(self):
        G = nx.Graph()
        u = (0.0, 0.0)
        n = (1.0, 1.0)
        v = (2.0, 0.0)

        G.add_node(u, tag="ref_vertex")
        G.add_node(n, tag="ref_vertex")
        G.add_node(v, tag="ref_vertex")
        G.add_edge(
            u,
            n,
            tag="ref_lines",
            geometry=LineString([u, n]),
            length=LineString([u, n]).length,
        )
        G.add_edge(
            n,
            v,
            tag="ref_lines",
            geometry=LineString([n, v]),
            length=LineString([n, v]).length,
        )

        bridge_with_straight_line(G, n)

        self.assertNotIn(n, G.nodes)
        self.assertTrue(G.has_edge(u, v))
        merged = G[u][v]["geometry"]
        self.assertEqual(list(merged.coords), [u, n, v])

    def test_prefer_vertices_ends_and_angles_removes_regular_real_vertex_near_end(self):
        G = nx.Graph()
        a = (0.0, 0.0)  # end
        b = (2.0, 0.0)  # regular interior vertex
        c = (4.0, 0.0)  # end

        G.add_node(a, tag="ref_vertex", is_reference_line_end=True)
        G.add_node(b, tag="ref_vertex")
        G.add_node(c, tag="ref_vertex", is_reference_line_end=True)
        G.add_edge(
            a,
            b,
            tag="ref_lines",
            geometry=LineString([a, b]),
            length=LineString([a, b]).length,
        )
        G.add_edge(
            b,
            c,
            tag="ref_lines",
            geometry=LineString([b, c]),
            length=LineString([b, c]).length,
        )

        cleaned = clean_pseudo_nodes_by_snap_strategy(
            G,
            snap_strategy=SnapStrategy.PREFER_VERTICES_ENDS_AND_ANGLES,
            distance_threshold=2.1,
            angle_threshold_degrees=150.0,
        )

        self.assertNotIn(b, cleaned.nodes)
        self.assertTrue(cleaned.has_edge(a, c))



    def test_pseudo_tag_interacts_with_snap_strategy_in_node_selection(self):
        G = nx.Graph()
        a = (0.0, 0.0)  # real end
        b = (2.0, 0.0)  # real interior
        c = (4.0, 0.0)  # real end
        p = (1.7, 0.0)  # pseudo node near start point

        G.add_node(a, tag="ref_vertex", is_reference_line_end=True)
        G.add_node(b, tag="ref_vertex")
        G.add_node(c, tag="ref_vertex", is_reference_line_end=True)
        G.add_node(p, tag="pseudo_ref_vertex")
        G.add_edge(a, b, geometry=LineString([a, b]), length=2.0, tag="ref_lines")
        G.add_edge(b, c, geometry=LineString([b, c]), length=2.0, tag="ref_lines")
        G.add_edge(p, b, geometry=LineString([p, b]), length=0.3, tag="interconnect")

        thematic = LineString([(1.6, 0.0), c])
        tolerance = 2.0

        result_no_pref = find_best_path_in_network(
            thematic,
            G,
            snap_strategy=SnapStrategy.NO_PREFERENCE,
            tolerance=tolerance,
        )
        result_ends_angles = find_best_path_in_network(
            thematic,
            G,
            snap_strategy=SnapStrategy.PREFER_VERTICES_ENDS_AND_ANGLES,
            tolerance=tolerance,
            angle_threshold_degrees=150.0,
        )

        self.assertIsNotNone(result_no_pref)
        self.assertIsNotNone(result_ends_angles)
        self.assertEqual(list(result_no_pref.coords), [p, b, c])
        self.assertEqual(list(result_ends_angles.coords), [a, b, c])
