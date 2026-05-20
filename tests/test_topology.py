import unittest
from unittest.mock import patch

import pytest
from shapely.geometry import LineString
from shapely.geometry import Polygon
from shapely import from_wkt

from brdr.aligner import Aligner
from brdr.be.grb.enums import GRBType
from brdr.be.grb.loader import GRBActualLoader
from brdr.configs import ProcessorConfig
from brdr.constants import REMARK_FIELD_NAME
from brdr.enums import ProcessRemark
from brdr.enums import SnapStrategy
from brdr.geometry_utils import safe_difference
from brdr.geometry_utils import safe_equals
from brdr.loader import DictLoader
from brdr.loader import GeoJsonLoader
from brdr.processor import NetworkGeometryProcessor
from brdr.processor import TopologyProcessor
from brdr.topo_utils import _dissolve_topo
from brdr.topo_utils import _generate_topo


class TestTopology(unittest.TestCase):
    def setUp(self):
        # Create a sample geometry for testing
        self.sample_aligner = Aligner()
        self.sample_geom = Polygon([(0, 0), (0, 10), (10, 10), (10, 0)])

    @pytest.mark.usefixtures("callback_grb_response")
    def test_topology(self):
        """
        Test if parameter preserve_topology is working"""
        # Initiate an Aligner
        geojson = {
            "type": "FeatureCollection",
            "name": "two",
            "crs": {
                "type": "name",
                "properties": {"name": "urn:ogc:def:crs:EPSG::31370"},
            },
            "features": [
                {
                    "type": "Feature",
                    "properties": {
                        "OIDN": 5763414,
                        "UIDN": 15089212,
                        "VERSIE": 1,
                        "BEGINDATUM": "2021-05-16",
                        "VERSDATUM": "2021-05-16",
                        "EINDDATUM": "2024-02-05",
                        "CAPAKEY": "24126B0049/00Z000",
                        "CANU": "49Z",
                        "FISCDATUM": "2022-01-01",
                        "BHRDR": 2,
                        "LBLBHRDR": "AAPD",
                        "NISCODE": "24062",
                        "BGNINV": 11,
                        "LBLBGNINV": "adpupdate",
                        "EINDINV": 11,
                        "LBLEINDINV": "adpupdate",
                        "BEWERK": 1,
                        "LBLBEWERK": "verwijderd",
                        "LENGTE": 117.62,
                        "OPPERVL": 398.32,
                    },
                    "geometry": {
                        "type": "MultiPolygon",
                        "coordinates": [
                            [
                                [
                                    [174193.075163975358009, 179501.127230677753687],
                                    [174193.069019973278046, 179501.100158676505089],
                                    [174183.54901996999979, 179504.310078680515289],
                                    [174184.360603965818882, 179506.715262681245804],
                                    [174192.349019974470139, 179530.390142697840929],
                                    [174197.229019977152348, 179552.460094712674618],
                                    [174199.52905198186636, 179551.690110713243484],
                                    [174204.018971979618073, 179550.170110709965229],
                                    [174196.94549997895956, 179518.199998687952757],
                                    [174196.899035975337029, 179517.990142688155174],
                                    [174193.623707972466946, 179503.546174678951502],
                                    [174193.075163975358009, 179501.127230677753687],
                                ]
                            ]
                        ],
                    },
                },
                {
                    "type": "Feature",
                    "properties": {
                        "OIDN": 3322239,
                        "UIDN": 15109808,
                        "VERSIE": 3,
                        "BEGINDATUM": "2011-12-08",
                        "VERSDATUM": "2021-05-16",
                        "EINDDATUM": "2024-02-06",
                        "CAPAKEY": "24126B0056/00F000",
                        "CANU": "56F",
                        "FISCDATUM": "2022-01-01",
                        "BHRDR": 2,
                        "LBLBHRDR": "AAPD",
                        "NISCODE": "24062",
                        "BGNINV": 11,
                        "LBLBGNINV": "adpupdate",
                        "EINDINV": 11,
                        "LBLEINDINV": "adpupdate",
                        "BEWERK": 3,
                        "LBLBEWERK": "geometriewijziging, niet beduidend",
                        "LENGTE": 290.35,
                        "OPPERVL": 2313.54,
                    },
                    "geometry": {
                        "type": "MultiPolygon",
                        "coordinates": [
                            [
                                [
                                    [174245.362332008779049, 179581.87654273211956],
                                    [174239.231516011059284, 179554.011070713400841],
                                    [174227.245403997600079, 179557.761854715645313],
                                    [174219.442075990140438, 179560.203838717192411],
                                    [174210.462427988648415, 179520.346942689269781],
                                    [174206.717275984585285, 179509.097086682915688],
                                    [174204.78287598490715, 179503.286270678043365],
                                    [174202.223963983356953, 179497.917822673916817],
                                    [174202.209243983030319, 179497.886910676956177],
                                    [174193.069019973278046, 179501.100158676505089],
                                    [174193.075163975358009, 179501.127230677753687],
                                    [174193.623707972466946, 179503.546174678951502],
                                    [174196.899035975337029, 179517.990142688155174],
                                    [174196.94549997895956, 179518.199998687952757],
                                    [174204.018971979618073, 179550.170110709965229],
                                    [174207.168987981975079, 179564.370046719908714],
                                    [174200.35887598246336, 179566.609854724258184],
                                    [174188.40898796916008, 179570.540030725300312],
                                    [174193.813019976019859, 179598.002366743981838],
                                    [174245.362332008779049, 179581.87654273211956],
                                ]
                            ]
                        ],
                    },
                },
            ],
        }

        processor = TopologyProcessor(config=ProcessorConfig(), feedback=None)
        aligner = Aligner(crs="EPSG:31370", processor=processor)
        loader = GeoJsonLoader(_input=geojson, id_property="CAPAKEY")
        aligner.load_thematic_data(loader)
        aligner.load_reference_data(
            GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
        )
        relevant_distance = 2
        aligner_result = aligner.process(relevant_distances=[relevant_distance])

        self.assertEqual(len(aligner_result.results), 2)
        dict_predictions_evaluated = aligner.evaluate()
        print(dict_predictions_evaluated)

    def test_topology_processor_duplicate_wkb_is_disambiguated_via_thematic_id(self):
        processor = TopologyProcessor(config=ProcessorConfig(), feedback=None)
        aligner = Aligner(crs="EPSG:31370", processor=processor)
        poly = Polygon([(0, 0), (0, 10), (10, 10), (10, 0)])
        thematic = {"a": poly, "b": poly}
        reference = {"r1": poly}
        aligner.load_thematic_data(DictLoader(thematic))
        aligner.load_reference_data(DictLoader(reference))

        result = aligner.process([1.0])
        assert len(result.results) == 2
        assert "a" in result.results
        assert "b" in result.results

    def test_topology_processor_duplicate_wkb_without_thematic_id_raises(self):
        processor = TopologyProcessor(config=ProcessorConfig(), feedback=None)
        aligner = Aligner(crs="EPSG:31370")
        poly = Polygon([(0, 0), (0, 10), (10, 10), (10, 0)])
        thematic = {"a": poly, "b": poly}
        reference = {"r1": poly}
        aligner.load_thematic_data(DictLoader(thematic))
        aligner.load_reference_data(DictLoader(reference))

        with self.assertRaises(ValueError):
            processor.process(
                input_geometry=poly,
                reference_data=aligner.reference_data,
                mitre_limit=2,
                correction_distance=0.01,
                relevant_distance=1.0,
                thematic_data=aligner.thematic_data,
            )

    def test_dissolve_topo_does_not_swallow_unexpected_runtime_error(self):
        aligner = Aligner(crs="EPSG:31370")
        poly = Polygon([(0, 0), (0, 10), (10, 10), (10, 0)])
        aligner.load_thematic_data(DictLoader({"a": poly}))
        thematic_geometries_to_process, topo = _generate_topo(aligner.thematic_data)
        arc_id = next(iter(thematic_geometries_to_process.keys()))
        process_results = {arc_id: {1.0: {"result": LineString([(0, 0), (10, 0)])}}}

        with patch(
            "brdr.topo_utils.longest_linestring_from_multilinestring",
            side_effect=RuntimeError("boom"),
        ):
            with self.assertRaises(RuntimeError):
                _dissolve_topo(
                    "a",
                    process_results,
                    poly,
                    thematic_geometries_to_process,
                    topo,
                    1.0,
                )

    def test_topology_arc_cache_invalidates_on_config_change(self):
        processor = TopologyProcessor(config=ProcessorConfig(), feedback=None)
        aligner = Aligner(crs="EPSG:31370")
        poly = Polygon([(0, 0), (0, 10), (10, 10), (10, 0)])
        aligner.load_thematic_data(DictLoader({"a": poly}))
        aligner.load_reference_data(DictLoader({"r1": poly}))

        with patch.object(
            NetworkGeometryProcessor,
            "process",
            autospec=True,
            wraps=NetworkGeometryProcessor.process,
        ) as mocked:
            # First run fills arc cache.
            processor.process(
                input_geometry=poly,
                thematic_id="a",
                reference_data=aligner.reference_data,
                mitre_limit=aligner.mitre_limit,
                correction_distance=aligner.correction_distance,
                relevant_distance=1.0,
                thematic_data=aligner.thematic_data,
            )
            calls_after_first = mocked.call_count
            assert calls_after_first > 0

            # Second identical run reuses cache.
            processor.process(
                input_geometry=poly,
                thematic_id="a",
                reference_data=aligner.reference_data,
                mitre_limit=aligner.mitre_limit,
                correction_distance=aligner.correction_distance,
                relevant_distance=1.0,
                thematic_data=aligner.thematic_data,
            )
            assert mocked.call_count == calls_after_first

            # Mutating behavior-driving config should invalidate cache key.
            processor.config.snap_strategy = SnapStrategy.NO_PREFERENCE
            processor.process(
                input_geometry=poly,
                thematic_id="a",
                reference_data=aligner.reference_data,
                mitre_limit=aligner.mitre_limit,
                correction_distance=aligner.correction_distance,
                relevant_distance=1.0,
                thematic_data=aligner.thematic_data,
            )
            assert mocked.call_count > calls_after_first

    def test_topology_processor_with_processing_area_no_overlap_returns_unchanged(self):
        processor = TopologyProcessor(config=ProcessorConfig(), feedback=None)
        aligner = Aligner(crs="EPSG:31370", processor=processor)
        thematic_geom = from_wkt("LINESTRING (0 0, 10 0)")
        thematic = {"t1": thematic_geom}
        reference = {"r1": from_wkt("LINESTRING (0 1, 10 1)")}
        scope = from_wkt("POLYGON ((20 20, 20 21, 21 21, 21 20, 20 20))")
        aligner.load_thematic_data(DictLoader(thematic))
        aligner.load_reference_data(DictLoader(reference))

        pr = aligner.process([1], processing_area=scope).results["t1"][1]
        assert safe_equals(pr["result"], thematic_geom)
        assert ProcessRemark.RESULT_UNCHANGED in pr["properties"].get(
            REMARK_FIELD_NAME, []
        )

    def test_topology_processor_with_processing_area_partial_only_changes_inside_scope(
        self,
    ):
        processor = TopologyProcessor(config=ProcessorConfig(), feedback=None)
        aligner = Aligner(crs="EPSG:31370", processor=processor)
        thematic_geom = from_wkt("POLYGON ((0 0, 0 6, 6 6, 6 0, 0 0))")
        thematic = {"t1": thematic_geom}
        reference = {"r1": from_wkt("POLYGON ((0 1, 0 7, 7 7, 7 1, 0 1))")}
        scope = from_wkt("POLYGON ((0 0, 0 6, 3 6, 3 0, 0 0))")
        aligner.load_thematic_data(DictLoader(thematic))
        aligner.load_reference_data(DictLoader(reference))

        result = aligner.process([1], processing_area=scope).results["t1"][1]["result"]
        outside = safe_difference(result, scope)
        expected_outside = safe_difference(thematic_geom, scope)

        assert not result.is_empty
        assert safe_equals(outside, expected_outside)

    def test_topology_processor_line_processing_area_partial_outside_contains_original_outside(
        self,
    ):
        processor = TopologyProcessor(config=ProcessorConfig(), feedback=None)
        aligner = Aligner(crs="EPSG:31370", processor=processor)
        thematic_geom = from_wkt("LINESTRING (0 0, 10 0)")
        thematic = {"t1": thematic_geom}
        reference = {"r1": from_wkt("LINESTRING (0 1, 10 1)")}
        scope = from_wkt("POLYGON ((0 -1, 0 2, 5 2, 5 -1, 0 -1))")
        aligner.load_thematic_data(DictLoader(thematic))
        aligner.load_reference_data(DictLoader(reference))

        result = aligner.process([1], processing_area=scope).results["t1"][1]["result"]
        outside = safe_difference(result, scope)
        expected_outside = safe_difference(thematic_geom, scope)

        assert not result.is_empty
        # For topology + lines, scoped processing may add extra out-of-scope segments,
        # but at minimum the original outside segment must be preserved.
        assert safe_difference(expected_outside, outside).is_empty
