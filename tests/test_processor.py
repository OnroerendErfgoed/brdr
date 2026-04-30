import unittest
from unittest.mock import patch

import numpy as np
from shapely import GeometryCollection
from shapely import from_wkt
from shapely.geometry import LineString
from shapely.geometry import Point

from brdr.aligner import Aligner
from brdr.configs import ProcessorConfig
from brdr.enums import OpenDomainStrategy
from brdr.geometry_utils import safe_difference
from brdr.graph_utils import _build_reference_segment_index
from brdr.loader import DictLoader
from brdr.processor import BaseProcessor
from brdr.processor import NetworkGeometryProcessor
from brdr.processor import AnchorGeometryProcessor


class _DummyProcessor(BaseProcessor):
    def process(self, **kwargs):
        raise NotImplementedError


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

    def test_postprocess_linear_od_exclude_ignores_non_polygon_reference(self):
        thematic = LineString([(0, 0), (4, 0)])
        preresult = LineString([(0, 0), (4, 0)])
        reference_union = LineString([(1, 0), (3, 0)])

        processor_as_is = _DummyProcessor(
            ProcessorConfig(od_strategy=OpenDomainStrategy.AS_IS)
        )
        out_as_is = processor_as_is._postprocess_preresult(
            geom_preresult=preresult,
            geom_thematic=thematic,
            relevant_intersection=GeometryCollection(),
            relevant_diff=GeometryCollection(),
            relevant_distance=0.5,
            reference_union=reference_union,
            mitre_limit=10,
            correction_distance=0.01,
        )

        processor_exclude = _DummyProcessor(
            ProcessorConfig(od_strategy=OpenDomainStrategy.EXCLUDE)
        )
        out_exclude = processor_exclude._postprocess_preresult(
            geom_preresult=preresult,
            geom_thematic=thematic,
            relevant_intersection=GeometryCollection(),
            relevant_diff=GeometryCollection(),
            relevant_distance=0.5,
            reference_union=reference_union,
            mitre_limit=10,
            correction_distance=0.01,
        )

        assert out_as_is["result"].length == 4.0
        assert out_exclude["result"].length == 4.0

    def test_postprocess_linear_od_exclude_clips_to_polygon_reference(self):
        thematic = LineString([(0, 0), (4, 0)])
        preresult = LineString([(0, 0), (4, 0)])
        reference_union = from_wkt("POLYGON ((1 -1, 3 -1, 3 1, 1 1, 1 -1))")

        processor_exclude = _DummyProcessor(
            ProcessorConfig(od_strategy=OpenDomainStrategy.EXCLUDE)
        )
        out_exclude = processor_exclude._postprocess_preresult(
            geom_preresult=preresult,
            geom_thematic=thematic,
            relevant_intersection=GeometryCollection(),
            relevant_diff=GeometryCollection(),
            relevant_distance=0.5,
            reference_union=reference_union,
            mitre_limit=10,
            correction_distance=0.01,
        )

        assert out_exclude["result"].length == 2.0

    def test_postprocess_point_od_exclude_ignores_non_polygon_reference(self):
        thematic = Point(0, 0)
        preresult = Point(0, 0)
        reference_union = LineString([(1, 0), (3, 0)])

        processor_as_is = _DummyProcessor(
            ProcessorConfig(od_strategy=OpenDomainStrategy.AS_IS)
        )
        out_as_is = processor_as_is._postprocess_preresult(
            geom_preresult=preresult,
            geom_thematic=thematic,
            relevant_intersection=GeometryCollection(),
            relevant_diff=GeometryCollection(),
            relevant_distance=0.5,
            reference_union=reference_union,
            mitre_limit=10,
            correction_distance=0.01,
        )

        processor_exclude = _DummyProcessor(
            ProcessorConfig(od_strategy=OpenDomainStrategy.EXCLUDE)
        )
        out_exclude = processor_exclude._postprocess_preresult(
            geom_preresult=preresult,
            geom_thematic=thematic,
            relevant_intersection=GeometryCollection(),
            relevant_diff=GeometryCollection(),
            relevant_distance=0.5,
            reference_union=reference_union,
            mitre_limit=10,
            correction_distance=0.01,
        )

        assert out_as_is["result"].geom_type == "Point"
        assert out_exclude["result"].geom_type == "Point"

    def test_postprocess_point_od_exclude_clips_to_polygon_reference(self):
        thematic = Point(0, 0)
        preresult = Point(0, 0)
        reference_union = from_wkt("POLYGON ((1 -1, 3 -1, 3 1, 1 1, 1 -1))")

        processor_exclude = _DummyProcessor(
            ProcessorConfig(od_strategy=OpenDomainStrategy.EXCLUDE)
        )
        out_exclude = processor_exclude._postprocess_preresult(
            geom_preresult=preresult,
            geom_thematic=thematic,
            relevant_intersection=GeometryCollection(),
            relevant_diff=GeometryCollection(),
            relevant_distance=0.5,
            reference_union=reference_union,
            mitre_limit=10,
            correction_distance=0.01,
        )

        assert out_exclude["result"].is_empty

    def test_network_od_exclude_vs_as_is_with_polygon_reference(self):
        thematic_dict = {"theme_id": from_wkt("LINESTRING (0 0, 10 0)")}
        reference_polygon = from_wkt("POLYGON ((2 -1, 8 -1, 8 1, 2 1, 2 -1))")
        reference_dict = {"ref_id": reference_polygon}

        aligner_exclude = Aligner(
            processor=NetworkGeometryProcessor(
                config=ProcessorConfig(od_strategy=OpenDomainStrategy.EXCLUDE)
            )
        )
        aligner_exclude.load_thematic_data(DictLoader(thematic_dict))
        aligner_exclude.load_reference_data(DictLoader(reference_dict))
        exclude_result = aligner_exclude.process([1]).results["theme_id"][1]["result"]

        aligner_as_is = Aligner(
            processor=NetworkGeometryProcessor(
                config=ProcessorConfig(od_strategy=OpenDomainStrategy.AS_IS)
            )
        )
        aligner_as_is.load_thematic_data(DictLoader(thematic_dict))
        aligner_as_is.load_reference_data(DictLoader(reference_dict))
        as_is_result = aligner_as_is.process([1]).results["theme_id"][1]["result"]

        assert safe_difference(exclude_result, reference_polygon).is_empty
        assert not safe_difference(as_is_result, reference_polygon).is_empty

    def test_network_all_od_strategies_with_polygon_reference(self):
        thematic_dict = {"theme_id": from_wkt("LINESTRING (0 0, 10 0)")}
        reference_dict = {"ref_id": from_wkt("POLYGON ((2 -1, 8 -1, 8 1, 2 1, 2 -1))")}

        for od_strategy in OpenDomainStrategy:
            aligner = Aligner(
                processor=NetworkGeometryProcessor(
                    config=ProcessorConfig(od_strategy=od_strategy)
                )
            )
            aligner.load_thematic_data(DictLoader(thematic_dict))
            aligner.load_reference_data(DictLoader(reference_dict))
            result = aligner.process([1]).results["theme_id"][1]["result"]
            assert result is not None

    def test_AnchorGeometryProcessor_aligns_linestring_to_reference_network(self):
        thematic_dict = {"theme_id": from_wkt("LINESTRING (0 0.2, 10 0.2)")}
        reference_dict = {"ref_id": from_wkt("LINESTRING (0 0, 10 0)")}

        processor = AnchorGeometryProcessor(config=ProcessorConfig())
        aligner = Aligner(processor=processor)
        aligner.load_thematic_data(DictLoader(thematic_dict))
        aligner.load_reference_data(DictLoader(reference_dict))

        result = aligner.process([1]).results["theme_id"][1]["result"]
        assert result is not None
        assert result.geom_type in {"LineString", "MultiLineString"}
        assert result.distance(reference_dict["ref_id"]) < 1e-8

    def test_anchorprocessor_uses_pathfinder_per_anchor_segment(self):
        processor = AnchorGeometryProcessor(
            config=ProcessorConfig(max_anchor_distance=5.0)
        )
        reference_lines = [
            LineString([(0, 0), (5, 0)]),
            LineString([(5, 0), (10, 0)]),
        ]
        graph = processor._build_reference_graph(reference_lines)
        line = LineString([(0, 0), (10, 0)])
        relevant_distance = 2.0

        def _echo_segment(segment_geom, graph_obj, **kwargs):
            assert isinstance(segment_geom, LineString)
            assert graph_obj is graph
            assert kwargs["snap_strategy"] == processor.config.snap_strategy
            assert kwargs["tolerance"] == relevant_distance
            return segment_geom

        with patch("brdr.processor.find_best_path_in_network", side_effect=_echo_segment) as mocked:
            result = processor._align_linestring_anchor_routing(
                line, graph, relevant_distance=relevant_distance
            )

        assert result is not None
        assert result.geom_type == "LineString"
        assert mocked.call_count == 2

    def test_anchorprocessor_anchor_guides_route_through_middle_corridor(self):
        processor = AnchorGeometryProcessor(config=ProcessorConfig())
        reference_lines = [
            LineString([(0, 0), (10, 0)]),
            LineString([(0, 10), (5, 10), (10, 10)]),
            LineString([(0, 0), (0, 10)]),
            LineString([(10, 0), (10, 10)]),
        ]
        graph = processor._build_reference_graph(reference_lines)
        input_line = LineString([(0, 0.2), (5, 9.8), (10, 0.2)])

        result = processor._align_linestring_anchor_routing(
            input_line, graph, relevant_distance=2.0
        )

        assert result is not None
        assert result.geom_type == "LineString"
        max_y = max(coord[1] for coord in result.coords)
        assert max_y >= 9.99

    def test_anchorprocessor_builds_directed_graph_with_oneway(self):
        processor = AnchorGeometryProcessor(
            config=ProcessorConfig(network_use_directed_graph=True)
        )
        ref_line = LineString([(0, 0), (10, 0)])
        reference_feature_records = [
            {
                "id": "ref1",
                "geometry": ref_line,
                "properties": {"oneway": "yes"},
            }
        ]
        direction_index = _build_reference_segment_index(
            reference_feature_records=reference_feature_records,
            oneway_field=processor.config.network_oneway_field,
            oneway_forward_values=processor.config.network_oneway_forward_values,
            oneway_reverse_values=processor.config.network_oneway_reverse_values,
        )

        graph = processor._build_reference_graph(
            [ref_line], precomputed_ref_direction_index=direction_index
        )

        assert graph.is_directed()
        assert graph.has_edge((0.0, 0.0), (10.0, 0.0))
        assert not graph.has_edge((10.0, 0.0), (0.0, 0.0))

    def test_anchorprocessor_routes_on_full_reference_after_local_anchor_match(self):
        thematic_dict = {"theme_id": from_wkt("LINESTRING (0 0.2, 10 0.2)")}
        reference_dict = {
            "left": from_wkt("LINESTRING (0 0, 0 5)"),
            "right": from_wkt("LINESTRING (10 0, 10 5)"),
            "top": from_wkt("LINESTRING (0 5, 10 5)"),
        }
        processor = AnchorGeometryProcessor(config=ProcessorConfig())
        aligner = Aligner(processor=processor)
        aligner.load_thematic_data(DictLoader(thematic_dict))
        aligner.load_reference_data(DictLoader(reference_dict))

        result = aligner.process([3]).results["theme_id"][3]["result"]
        assert result is not None
        assert result.geom_type in {"LineString", "MultiLineString"}
        top_line = from_wkt("LINESTRING (0 5, 10 5)")
        overlap_with_top = result.intersection(top_line)
        assert not overlap_with_top.is_empty

    def test_anchorprocessor_keeps_unmatched_anchor_point_in_route(self):
        processor = AnchorGeometryProcessor(config=ProcessorConfig())
        reference_lines = [LineString([(0, 0), (10, 0)])]
        graph = processor._build_reference_graph(reference_lines)
        # Mid anchor at (5, 2) has no candidate within relevant_distance=1 and must stay.
        input_line = LineString([(0, 0), (5, 2), (10, 0)])

        result = processor._align_linestring_anchor_routing(
            input_line,
            graph,
            relevant_distance=1.0,
            anchor_selection_graph=graph,
        )

        assert result is not None
        assert result.geom_type == "LineString"
        assert (5.0, 2.0) in list(result.coords)

    def test_anchorprocessor_anchor_points_include_sharp_vertex(self):
        processor = AnchorGeometryProcessor(config=ProcessorConfig())
        line = LineString([(0, 0), (5, 5), (10, 0)])

        anchors = processor._anchor_points(line, relevant_distance=2.0)
        anchor_coords = {(round(p.x, 6), round(p.y, 6)) for p in anchors}

        assert (0.0, 0.0) in anchor_coords
        assert (5.0, 5.0) in anchor_coords
        assert (10.0, 0.0) in anchor_coords

    def test_anchorprocessor_anchor_points_include_periodic_points_on_long_line(self):
        processor = AnchorGeometryProcessor(config=ProcessorConfig())
        line = LineString([(0, 0), (20, 0)])

        anchors = processor._anchor_points(line, relevant_distance=2.0)
        anchor_coords = {(round(p.x, 6), round(p.y, 6)) for p in anchors}
        # start + periodic at 10 + end
        assert (0.0, 0.0) in anchor_coords
        assert (10.0, 0.0) in anchor_coords
        assert (20.0, 0.0) in anchor_coords

    def test_anchorprocessor_anchor_points_use_configured_max_anchor_distance(self):
        processor = AnchorGeometryProcessor(
            config=ProcessorConfig(max_anchor_distance=6.0)
        )
        line = LineString([(0, 0), (20, 0)])

        anchors = processor._anchor_points(line, relevant_distance=2.0)
        anchor_coords = {(round(p.x, 6), round(p.y, 6)) for p in anchors}
        assert (6.0, 0.0) in anchor_coords
        assert (12.0, 0.0) in anchor_coords
        assert (18.0, 0.0) in anchor_coords

