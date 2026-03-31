import unittest

import numpy as np
from shapely import GeometryCollection
from shapely import from_wkt
from shapely.geometry import Point
from shapely.geometry import LineString

from brdr.aligner import Aligner
from brdr.configs import ProcessorConfig
from brdr.enums import OpenDomainStrategy
from brdr.geometry_utils import safe_difference
from brdr.loader import DictLoader
from brdr.processor import BaseProcessor
from brdr.processor import NetworkGeometryProcessor


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
        reference_union = from_wkt(
            "POLYGON ((1 -1, 3 -1, 3 1, 1 1, 1 -1))"
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
        reference_union = from_wkt(
            "POLYGON ((1 -1, 3 -1, 3 1, 1 1, 1 -1))"
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
        reference_dict = {
            "ref_id": from_wkt("POLYGON ((2 -1, 8 -1, 8 1, 2 1, 2 -1))")
        }

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
