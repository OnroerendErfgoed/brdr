import unittest

import networkx as nx
from shapely.geometry import LineString, MultiLineString

from brdr.configs import ProcessorConfig
from brdr.enums import SnapStrategy
from brdr.graph_utils import build_custom_network, find_best_path_in_network
from brdr.processor import DirectedNetworkGeometryProcessor


class TestDirectedNetwork(unittest.TestCase):
    def _build_graph(
        self,
        thematic_line,
        reference_line,
        *,
        directed,
        oneway_value=None,
        allow_connectors=False,
        penalty=10.0,
    ):
        records = [
            {
                "id": "ref1",
                "geometry": reference_line,
                "properties": {} if oneway_value is None else {"oneway": oneway_value},
            }
        ]
        return build_custom_network(
            input_geometry=thematic_line,
            theme_multiline=MultiLineString([thematic_line]),
            reference=reference_line,
            reference_intersection=reference_line,
            reference_feature_records=records,
            relevant_distance=1.0,
            snap_dist=0.001,
            directed=directed,
            oneway_field="oneway",
            oneway_forward_values=("yes", "1", "true", "forward"),
            oneway_reverse_values=("-1", "reverse", "backward"),
            allow_connector_edges_when_directed=allow_connectors,
            directed_connector_penalty_factor=penalty,
        )

    def test_undirected_mode_keeps_reverse_path_possible(self):
        thematic = LineString([(10.0, 0.0), (0.0, 0.0)])
        reference = LineString([(0.0, 0.0), (10.0, 0.0)])
        graph = self._build_graph(
            thematic, reference, directed=False, oneway_value="yes"
        )
        result = find_best_path_in_network(
            thematic,
            graph,
            snap_strategy=SnapStrategy.NO_PREFERENCE,
            tolerance=1.0,
        )
        self.assertIsNotNone(result)

    def test_directed_mode_blocks_reverse_on_forward_oneway(self):
        thematic = LineString([(10.0, 0.0), (0.0, 0.0)])
        reference = LineString([(0.0, 0.0), (10.0, 0.0)])
        graph = self._build_graph(thematic, reference, directed=True, oneway_value="yes")
        result = find_best_path_in_network(
            thematic,
            graph,
            snap_strategy=SnapStrategy.NO_PREFERENCE,
            tolerance=1.0,
        )
        self.assertIsNone(result)

    def test_directed_mode_allows_forward_on_forward_oneway(self):
        thematic = LineString([(0.0, 0.0), (10.0, 0.0)])
        reference = LineString([(0.0, 0.0), (10.0, 0.0)])
        graph = self._build_graph(thematic, reference, directed=True, oneway_value="yes")
        result = find_best_path_in_network(
            thematic,
            graph,
            snap_strategy=SnapStrategy.NO_PREFERENCE,
            tolerance=1.0,
        )
        self.assertIsNotNone(result)

    def test_directed_mode_reverse_value_allows_reverse(self):
        thematic = LineString([(10.0, 0.0), (0.0, 0.0)])
        reference = LineString([(0.0, 0.0), (10.0, 0.0)])
        graph = self._build_graph(thematic, reference, directed=True, oneway_value="-1")
        result = find_best_path_in_network(
            thematic,
            graph,
            snap_strategy=SnapStrategy.NO_PREFERENCE,
            tolerance=1.0,
        )
        self.assertIsNotNone(result)

    def test_closed_ring_in_directed_mode_uses_undirected_cycle_search(self):
        graph = nx.DiGraph()
        cycle = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
        for i in range(len(cycle)):
            u = cycle[i]
            v = cycle[(i + 1) % len(cycle)]
            geom = LineString([u, v])
            graph.add_edge(u, v, geometry=geom, length=geom.length, tag="ref_lines")

        thematic_ring = LineString(cycle + [cycle[0]])
        result = find_best_path_in_network(
            thematic_ring,
            graph,
            snap_strategy=SnapStrategy.NO_PREFERENCE,
            tolerance=1.0,
        )
        self.assertIsNotNone(result)
        self.assertTrue(result.is_ring)

    def test_directed_connector_edges_are_penalized(self):
        thematic = LineString([(0.0, 0.0), (1.0, 0.0)])
        reference = LineString([(2.0, 0.0), (3.0, 0.0)])
        penalty = 7.0
        graph = self._build_graph(
            thematic,
            reference,
            directed=True,
            oneway_value=None,
            allow_connectors=True,
            penalty=penalty,
        )
        interconnect_lengths = [
            d.get("length")
            for _, _, d in graph.edges(data=True)
            if d.get("tag") == "interconnect_lines"
        ]
        self.assertTrue(interconnect_lengths)
        self.assertTrue(all(length >= penalty for length in interconnect_lengths))

    def test_directed_processor_forces_directed_without_mutating_input_config(self):
        base_config = ProcessorConfig(network_use_directed_graph=False)
        processor = DirectedNetworkGeometryProcessor(config=base_config)
        self.assertFalse(base_config.network_use_directed_graph)
        self.assertTrue(processor.config.network_use_directed_graph)

    def test_directed_processor_blocks_reverse_on_forward_oneway(self):
        processor = DirectedNetworkGeometryProcessor(
            config=ProcessorConfig(network_use_directed_graph=False)
        )
        thematic = LineString([(10.0, 0.0), (0.0, 0.0)])
        reference = LineString([(0.0, 0.0), (10.0, 0.0)])
        records = [{"id": "ref1", "geometry": reference, "properties": {"oneway": "yes"}}]
        result = processor._get_processed_network_path(
            input_geometry=thematic,
            reference=reference,
            reference_feature_records=records,
            reference_intersection=reference,
            thematic_difference=MultiLineString([thematic]),
            relevant_distance=1.0,
            correction_distance=0.001,
        )
        self.assertIsNone(result)
