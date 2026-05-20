import unittest
from unittest.mock import patch

import networkx as nx
import numpy as np
from shapely import GeometryCollection, Polygon
from shapely import from_wkt
from shapely.geometry import LineString
from shapely.geometry import Point

from brdr.aligner import Aligner
from brdr.configs import ProcessorConfig
from brdr.constants import REMARK_FIELD_NAME
from brdr.enums import OpenDomainStrategy
from brdr.enums import ProcessRemark
from brdr.geometry_utils import safe_difference
from brdr.graph_utils import _build_reference_segment_index
from brdr.graph_utils import find_best_circle_path
from brdr.loader import DictLoader
from brdr.processor import BaseProcessor, DirectedNetworkGeometryProcessor
from brdr.processor import DieussaertGeometryProcessor
from brdr.processor import NetworkGeometryProcessor
from brdr.processor import AnchorGeometryProcessor
from brdr.processor import DirectedAnchorGeometryProcessor
from brdr.processor import AlignerGeometryProcessor


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

    def test_process_with_processing_area_helper_no_overlap_returns_unchanged(self):
        processor = _DummyProcessor(ProcessorConfig())
        input_geometry = LineString([(0, 0), (10, 0)])
        processing_area = from_wkt("POLYGON ((20 20, 20 21, 21 21, 21 20, 20 20))")

        result = processor._process_with_processing_area(
            input_geometry=input_geometry,
            processing_area=processing_area,
            reference_union=GeometryCollection(),
            relevant_distance=1.0,
            mitre_limit=10.0,
            correction_distance=0.01,
            scoped_processor=lambda g: {"result": g, "properties": {REMARK_FIELD_NAME: []}},
        )

        assert result["result"].equals(input_geometry)
        assert ProcessRemark.RESULT_UNCHANGED in result["properties"][REMARK_FIELD_NAME]

    def test_process_with_processing_area_helper_runs_postprocess_after_merge(self):
        processor = _DummyProcessor(ProcessorConfig())
        input_geometry = from_wkt("POLYGON ((0 0, 0 4, 10 4, 10 0, 0 0))")
        processing_area = from_wkt("POLYGON ((0 0, 0 4, 5 4, 5 0, 0 0))")

        scoped_result = from_wkt("POLYGON ((0 0, 0 4, 4.8 4, 4.8 0, 0 0))")
        with patch.object(
            processor,
            "_postprocess_preresult",
            wraps=processor._postprocess_preresult,
        ) as postprocess_spy:
            result = processor._process_with_processing_area(
                input_geometry=input_geometry,
                processing_area=processing_area,
                reference_union=GeometryCollection(),
                relevant_distance=1.0,
                mitre_limit=10.0,
                correction_distance=0.01,
                scoped_processor=lambda _g: {
                    "result": scoped_result,
                    "properties": {REMARK_FIELD_NAME: []},
                },
            )

        assert postprocess_spy.call_count == 1
        assert result["result"] is not None
        assert not result["result"].is_empty

    def test_aligner_processor_with_processing_area_forces_network_for_polygons(self):
        thematic = {"t1": from_wkt("POLYGON ((0 0, 0 6, 6 6, 6 0, 0 0))")}
        reference = {"r1": from_wkt("POLYGON ((0 1, 0 7, 7 7, 7 1, 0 1))")}
        scope = from_wkt("POLYGON ((0 0, 0 6, 3 6, 3 0, 0 0))")
        aligner = Aligner(processor=AlignerGeometryProcessor(config=ProcessorConfig()))
        aligner.load_thematic_data(DictLoader(thematic))
        aligner.load_reference_data(DictLoader(reference))

        with patch.object(
            DieussaertGeometryProcessor,
            "process",
            autospec=True,
            wraps=DieussaertGeometryProcessor.process,
        ) as dieussaert_process_spy:
            result = aligner.process([1], processing_area=scope).results["t1"][1]["result"]

        assert result is not None
        assert not result.is_empty
        assert dieussaert_process_spy.call_count == 0

    def test_aligner_processor_polygon_processing_area_uses_generic_scoped_flow(self):
        thematic = {"t1": from_wkt("POLYGON ((0 0, 0 6, 6 6, 6 0, 0 0))")}
        reference = {"r1": from_wkt("POLYGON ((0 1, 0 7, 7 7, 7 1, 0 1))")}
        scope = from_wkt("POLYGON ((0 0, 0 6, 3 6, 3 0, 0 0))")
        aligner = Aligner(processor=AlignerGeometryProcessor(config=ProcessorConfig()))
        aligner.load_thematic_data(DictLoader(thematic))
        aligner.load_reference_data(DictLoader(reference))

        result = aligner.process([1], processing_area=scope).results["t1"][1]["result"]
        assert result is not None
        assert not result.is_empty

    def test_aligner_processor_non_polygon_processing_area_uses_generic_scoped_flow(self):
        thematic = {"t1": from_wkt("POLYGON ((0 0, 0 6, 6 6, 6 0, 0 0))")}
        reference = {"r1": from_wkt("POLYGON ((0 1, 0 7, 7 7, 7 1, 0 1))")}
        # Non-polygon scope: generic scoped flow should handle this.
        scope = from_wkt("LINESTRING (0 0, 6 6)")
        aligner = Aligner(processor=AlignerGeometryProcessor(config=ProcessorConfig()))
        aligner.load_thematic_data(DictLoader(thematic))
        aligner.load_reference_data(DictLoader(reference))

        result = aligner.process([1], processing_area=scope).results["t1"][1]["result"]
        assert result is not None
        assert not result.is_empty

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

    def test_postprocess_invalid_preresult_none_returns_empty(self):
        thematic = from_wkt("POLYGON ((0 0, 0 2, 2 2, 2 0, 0 0))")
        processor = _DummyProcessor(ProcessorConfig())
        out = processor._postprocess_preresult(
            geom_preresult=None,
            geom_thematic=thematic,
            relevant_intersection=GeometryCollection(),
            relevant_diff=GeometryCollection(),
            relevant_distance=1.0,
            reference_union=GeometryCollection(),
            mitre_limit=10,
            correction_distance=0.01,
        )
        assert out["result"].is_empty
        assert ProcessRemark.INVALID_EMPTY_RETURNED in out["properties"].get(
            REMARK_FIELD_NAME, []
        )

    def test_postprocess_invalid_preresult_type_mismatch_returns_empty(self):
        thematic = from_wkt("POLYGON ((0 0, 0 2, 2 2, 2 0, 0 0))")
        processor = _DummyProcessor(ProcessorConfig())
        out = processor._postprocess_preresult(
            geom_preresult=LineString([(0, 0), (1, 0)]),
            geom_thematic=thematic,
            relevant_intersection=GeometryCollection(),
            relevant_diff=GeometryCollection(),
            relevant_distance=1.0,
            reference_union=GeometryCollection(),
            mitre_limit=10,
            correction_distance=0.01,
        )
        assert out["result"].is_empty
        assert ProcessRemark.INVALID_EMPTY_RETURNED in out["properties"].get(
            REMARK_FIELD_NAME, []
        )

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

    def test_dieussaert_calculate_geom_exclusion_returns_empty_polygons(self):
        processor = DieussaertGeometryProcessor(
            config=ProcessorConfig(threshold_exclusion_percentage=5.0)
        )
        geom_intersection = from_wkt("POLYGON ((0 0, 0 1, 1 1, 1 0, 0 0))")
        geom_reference = from_wkt("POLYGON ((0 0, 0 10, 10 10, 10 0, 0 0))")
        input_inner = from_wkt("POLYGON EMPTY")

        geom, rel_i, rel_d = processor._calculate_geom_by_intersection_and_reference(
            geom_intersection=geom_intersection,
            geom_reference=geom_reference,
            input_geometry_inner=input_inner,
            is_open_domain=False,
            buffer_distance=0.5,
            mitre_limit=5.0,
        )

        assert geom.geom_type == "Polygon" and geom.is_empty
        assert rel_i.geom_type == "Polygon" and rel_i.is_empty
        assert rel_d.geom_type == "Polygon" and rel_d.is_empty

    def test_dieussaert_calculate_geom_full_inclusion_returns_reference(self):
        processor = DieussaertGeometryProcessor(config=ProcessorConfig())
        geom_reference = from_wkt("POLYGON ((0 0, 0 4, 4 4, 4 0, 0 0))")
        input_inner = from_wkt("POLYGON EMPTY")

        geom, rel_i, rel_d = processor._calculate_geom_by_intersection_and_reference(
            geom_intersection=geom_reference,
            geom_reference=geom_reference,
            input_geometry_inner=input_inner,
            is_open_domain=False,
            buffer_distance=0.5,
            mitre_limit=5.0,
        )

        assert geom.equals(geom_reference)
        assert rel_i.equals(geom_reference)
        assert rel_d.geom_type == "Polygon" and rel_d.is_empty

    def test_dieussaert_threshold_overlap_minus_one_returns_intersection(self):
        processor = DieussaertGeometryProcessor(
            config=ProcessorConfig(threshold_overlap_percentage=-1)
        )
        geom_reference = from_wkt("POLYGON ((0 0, 0 10, 10 10, 10 0, 0 0))")
        geom_intersection = from_wkt("POLYGON ((2 2, 2 8, 8 8, 8 2, 2 2))")
        input_inner = from_wkt("POLYGON EMPTY")

        with patch("brdr.processor.buffer_neg", return_value=Polygon()):
            geom, _, _ = processor._calculate_geom_by_intersection_and_reference(
                geom_intersection=geom_intersection,
                geom_reference=geom_reference,
                input_geometry_inner=input_inner,
                is_open_domain=False,
                buffer_distance=0.5,
                mitre_limit=5.0,
            )

        assert geom.equals(geom_intersection)

    def test_dieussaert_open_domain_branch_uses_snap_geometry(self):
        processor = DieussaertGeometryProcessor(config=ProcessorConfig())
        geom_intersection = from_wkt("POLYGON ((2 2, 2 8, 8 8, 8 2, 2 2))")
        input_inner = from_wkt("POLYGON EMPTY")

        with patch("brdr.processor.buffer_neg", return_value=Polygon()):
            with patch(
                "brdr.processor.snap_geometry_to_reference",
                side_effect=lambda g, *_args, **_kwargs: g,
            ) as mocked_snap:
                geom, _, _ = processor._calculate_geom_by_intersection_and_reference(
                    geom_intersection=geom_intersection,
                    geom_reference=GeometryCollection(),  # area == 0 => OD reference
                    input_geometry_inner=input_inner,
                    is_open_domain=True,
                    buffer_distance=0.5,
                    mitre_limit=5.0,
                )

        assert mocked_snap.called
        assert geom.equals(geom_intersection)

    def test_dieussaert_partial_snapping_flag_controls_snap_call(self):
        geom_reference = from_wkt("POLYGON ((0 0, 0 10, 10 10, 10 0, 0 0))")
        geom_intersection = from_wkt("POLYGON ((4 4, 4 6, 6 6, 6 4, 4 4))")
        input_inner = from_wkt("POLYGON EMPTY")

        processor_no_snap = DieussaertGeometryProcessor(
            config=ProcessorConfig(partial_snapping=False)
        )
        with patch(
            "brdr.processor.snap_geometry_to_reference",
            side_effect=lambda g, *_args, **_kwargs: g,
        ) as mocked_snap_false:
            processor_no_snap._calculate_geom_by_intersection_and_reference(
                geom_intersection=geom_intersection,
                geom_reference=geom_reference,
                input_geometry_inner=input_inner,
                is_open_domain=False,
                buffer_distance=0.5,
                mitre_limit=5.0,
            )
        calls_without = mocked_snap_false.call_count

        processor_with_snap = DieussaertGeometryProcessor(
            config=ProcessorConfig(partial_snapping=True)
        )
        with patch(
            "brdr.processor.snap_geometry_to_reference",
            side_effect=lambda g, *_args, **_kwargs: g,
        ) as mocked_snap_true:
            processor_with_snap._calculate_geom_by_intersection_and_reference(
                geom_intersection=geom_intersection,
                geom_reference=geom_reference,
                input_geometry_inner=input_inner,
                is_open_domain=False,
                buffer_distance=0.5,
                mitre_limit=5.0,
            )
        calls_with = mocked_snap_true.call_count

        assert calls_with >= calls_without
        assert calls_with > 0

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

    def test_directed_anchor_processor_forces_directed_without_mutating_input_config(self):
        base_config = ProcessorConfig(network_use_directed_graph=False)
        processor = DirectedAnchorGeometryProcessor(config=base_config)
        assert not base_config.network_use_directed_graph
        assert processor.config.network_use_directed_graph

    def test_directed_anchor_no_path_returns_original_line_instead_of_direct_jump(self):
        processor = DirectedAnchorGeometryProcessor(config=ProcessorConfig())
        # One-way reference edge: only (0,0) -> (10,0)
        ref_line = LineString([(0, 0), (10, 0)])
        direction_index = _build_reference_segment_index(
            reference_feature_records=[
                {
                    "id": "ref1",
                    "geometry": ref_line,
                    "properties": {"oneway": "yes"},
                }
            ],
            oneway_field=processor.config.network_oneway_field,
            oneway_forward_values=processor.config.network_oneway_forward_values,
            oneway_reverse_values=processor.config.network_oneway_reverse_values,
        )
        graph = processor._build_reference_graph(
            [ref_line], precomputed_ref_direction_index=direction_index
        )
        input_line = LineString([(10, 0), (0, 0)])

        result = processor._align_linestring_anchor_routing(
            input_line,
            graph,
            relevant_distance=1.0,
            anchor_selection_graph=graph,
        )

        assert result is not None
        assert result.equals(input_line)

    def test_directed_anchor_builds_direction_index_from_all_reference_candidates(self):
        thematic_dict = {"theme_id": from_wkt("LINESTRING (0 0, 10 0)")}
        reference_dict = {
            "near": from_wkt("LINESTRING (0 0, 10 0)"),
            "far": from_wkt("LINESTRING (100 100, 110 100)"),
        }
        processor = DirectedAnchorGeometryProcessor(config=ProcessorConfig())
        aligner = Aligner(processor=processor)
        aligner.load_thematic_data(DictLoader(thematic_dict))
        aligner.load_reference_data(DictLoader(reference_dict))

        seen_ids = []
        original_builder = __import__(
            "brdr.processor", fromlist=["_build_reference_segment_index"]
        )._build_reference_segment_index

        def _capture_index(*, reference_feature_records, **kwargs):
            seen_ids.append({rec["id"] for rec in reference_feature_records})
            return original_builder(
                reference_feature_records=reference_feature_records, **kwargs
            )

        with patch("brdr.processor._build_reference_segment_index", side_effect=_capture_index):
            aligner.processor.process(
                input_geometry=thematic_dict["theme_id"],
                reference_data=aligner.reference_data,
                mitre_limit=aligner.mitre_limit,
                correction_distance=aligner.correction_distance,
                relevant_distance=1.0,
                reference_candidates=["near", "far"],
            )

        assert seen_ids
        assert {"near", "far"} in seen_ids

    def test_directed_network_no_path_returns_original_line_without_component_fallback(self):
        processor = DirectedNetworkGeometryProcessor(config=ProcessorConfig())
        thematic = LineString([(10, 0), (0, 0)])
        reference_dict = {
            "ref1": from_wkt("LINESTRING (0 0, 10 0)"),
        }
        aligner = Aligner(processor=processor)
        aligner.load_thematic_data(DictLoader({"theme_id": thematic}))
        aligner.load_reference_data(DictLoader(reference_dict))

        result = aligner.process([1.0]).results["theme_id"][1.0]["result"]
        assert result is not None
        # With strict directed routing and no reverse path, behavior should no longer
        # create a synthetic reverse connector route through fallback component linking.
        assert result.equals(thematic)

    def test_directed_network_builds_direction_index_from_all_reference_candidates(self):
        processor = DirectedNetworkGeometryProcessor(config=ProcessorConfig())
        thematic_dict = {"theme_id": from_wkt("LINESTRING (0 0, 10 0)")}
        reference_dict = {
            "near": from_wkt("LINESTRING (0 0, 10 0)"),
            "far": from_wkt("LINESTRING (100 100, 110 100)"),
        }
        aligner = Aligner(processor=processor)
        aligner.load_thematic_data(DictLoader(thematic_dict))
        aligner.load_reference_data(DictLoader(reference_dict))

        seen_ids = []
        original_builder = __import__(
            "brdr.processor", fromlist=["_build_reference_segment_index"]
        )._build_reference_segment_index

        def _capture_index(*, reference_feature_records, **kwargs):
            seen_ids.append({rec["id"] for rec in reference_feature_records})
            return original_builder(
                reference_feature_records=reference_feature_records, **kwargs
            )

        with patch("brdr.processor._build_reference_segment_index", side_effect=_capture_index):
            aligner.processor.process(
                input_geometry=thematic_dict["theme_id"],
                reference_data=aligner.reference_data,
                mitre_limit=aligner.mitre_limit,
                correction_distance=aligner.correction_distance,
                relevant_distance=1.0,
                reference_candidates=["near", "far"],
            )

        assert seen_ids
        assert {"near", "far"} in seen_ids

    def test_build_custom_network_unconnected_distance_scales_with_relevant_distance(self):
        from brdr.graph_utils import build_custom_network

        input_geometry = LineString([(0, 0), (10, 0)])
        theme_multiline = LineString([(0, 0), (10, 0)])
        reference = LineString([(0, 0), (10, 0)])
        reference_intersection = LineString([(0, 0), (10, 0)])

        max_seen = {"value": None}

        def _capture_connect_unconnected_greedy(G, max_spatial_dist=50, detour_ratio=3.0):
            max_seen["value"] = max_spatial_dist
            return G

        with patch(
            "brdr.graph_utils.connect_unconnected_greedy",
            side_effect=_capture_connect_unconnected_greedy,
        ):
            build_custom_network(
                input_geometry=input_geometry,
                theme_multiline=theme_multiline,
                reference=reference,
                reference_intersection=reference_intersection,
                relevant_distance=2.0,
                directed=False,
            )

        assert max_seen["value"] == 50

    def test_find_best_circle_path_single_polygon_low_overlap_returns_original_ring(self):
        graph = nx.Graph()
        # Candidate polygon around x=100 (far from original around x=0..10)
        coords = [(100, 0), (100, 10), (110, 10), (110, 0), (100, 0)]
        for i in range(len(coords) - 1):
            a = coords[i]
            b = coords[i + 1]
            graph.add_edge(a, b, geometry=LineString([a, b]), length=1.0)

        original_ring = LineString([(0, 0), (0, 10), (10, 10), (10, 0), (0, 0)])
        out = find_best_circle_path(graph, original_ring)
        assert out is original_ring

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
        processor = AnchorGeometryProcessor(
            config=ProcessorConfig(angle_threshold_degrees=180.0)
        )
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

        anchors = processor._anchor_points(line)
        anchor_coords = {(round(p.x, 6), round(p.y, 6)) for p in anchors}

        assert (0.0, 0.0) in anchor_coords
        assert (5.0, 5.0) in anchor_coords
        assert (10.0, 0.0) in anchor_coords

    def test_anchorprocessor_anchor_points_include_periodic_points_on_long_line(self):
        processor = AnchorGeometryProcessor(
            config=ProcessorConfig(max_anchor_distance=10.0)
        )
        line = LineString([(0, 0), (20, 0)])

        anchors = processor._anchor_points(line)
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

        anchors = processor._anchor_points(line)
        anchor_coords = {(round(p.x, 6), round(p.y, 6)) for p in anchors}
        assert (6.0, 0.0) in anchor_coords
        assert (12.0, 0.0) in anchor_coords
        assert (18.0, 0.0) in anchor_coords
