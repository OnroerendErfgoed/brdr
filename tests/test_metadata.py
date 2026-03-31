import json
import unittest

from shapely import from_wkt

from brdr.aligner import (
    Aligner,
    _get_metadata_observations_from_process_result,
    _reverse_metadata_observations_to_brdr_observation,
)
from brdr.configs import AlignerConfig
from brdr.loader import DictLoader


def test_metadata():
    config = AlignerConfig(
        log_metadata=True,
        add_observations=True,
    )
    aligner = Aligner(config=config)
    aligner.load_thematic_data(
        DictLoader({"theme_id_1": from_wkt("POLYGON ((0 0, 0 9, 5 10, 10 0, 0 0))")})
    )
    aligner.load_reference_data(
        DictLoader({"ref_id_1": from_wkt("POLYGON ((0 1, 0 10,8 10,10 1,0 1))")})
    )
    result = aligner.process([1])
    assert True  # TODO


class TestReverseMetadata(unittest.TestCase):

    def setUp(self):
        """Mock data to re-use"""
        metadata_string = '{"actuation": {"id": "urn:uuid:ea9f6199-9852-4dce-96b0-48428d928e9e", "type": "sosa:Actuation", "reference_geometries": [{"id": "urn:uuid:6b46bf10-1c72-42d3-8c91-175d63519caa", "type": "geo:Polygon", "version_date": "2026-02-04", "derived_from": {"id": "24126B0030/00G000", "type": "geo:Feature", "source": "https://geo.api.vlaanderen.be/GRB/ogc/features/collections/ADP"}}, {"id": "urn:uuid:7fed607b-7215-409a-802e-d55228831a07", "type": "geo:Polygon", "version_date": "2026-02-04", "derived_from": {"id": "24126B0030/00R000", "type": "geo:Feature", "source": "https://geo.api.vlaanderen.be/GRB/ogc/features/collections/ADP"}}, {"id": "urn:uuid:62690bd4-8e1e-4101-a3bf-27a88eb61518", "type": "geo:Polygon", "version_date": "2026-02-04", "derived_from": {"id": "24126B0029/00T000", "type": "geo:Feature", "source": "https://geo.api.vlaanderen.be/GRB/ogc/features/collections/ADP"}}, {"id": "urn:uuid:643b9345-815b-4e0f-bab8-9aa500811705", "type": "geo:Polygon", "version_date": "2026-02-04", "derived_from": {"id": "24126B0029/00C002", "type": "geo:Feature", "source": "https://geo.api.vlaanderen.be/GRB/ogc/features/collections/ADP"}}, {"id": "urn:uuid:59533634-a68c-4a1f-be35-c8bdb05db3ad", "type": "geo:Polygon", "version_date": "2026-02-04", "derived_from": {"id": "24126B0029/00V000", "type": "geo:Feature", "source": "https://geo.api.vlaanderen.be/GRB/ogc/features/collections/ADP"}}, {"id": "urn:uuid:863c0d38-4a8b-4184-a171-d6b0cd4ddc13", "type": "geo:Polygon", "version_date": "2026-02-04", "derived_from": {"id": "24126B0029/00S000", "type": "geo:Feature", "source": "https://geo.api.vlaanderen.be/GRB/ogc/features/collections/ADP"}}, {"id": "urn:uuid:302ba351-8ad2-4139-832b-03e4445dd012", "type": "geo:Polygon", "version_date": "2026-02-04", "derived_from": {"id": "24126B0029/00R000", "type": "geo:Feature", "source": "https://geo.api.vlaanderen.be/GRB/ogc/features/collections/ADP"}}, {"id": "urn:uuid:0b0ab232-f78c-4378-8e1c-1be1b36b6718", "type": "geo:Polygon", "version_date": "2026-02-04", "derived_from": {"id": "24126B0029/00W000", "type": "geo:Feature", "source": "https://geo.api.vlaanderen.be/GRB/ogc/features/collections/ADP"}}, {"id": "urn:uuid:13b4c5f3-4755-4cb0-8e8c-08e7fe25ccf5", "type": "geo:Polygon", "version_date": "2026-02-04", "derived_from": {"id": "24126B0027/00X002", "type": "geo:Feature", "source": "https://geo.api.vlaanderen.be/GRB/ogc/features/collections/ADP"}}], "changes": "geo:hasGeometry", "sosa:hasFeatureOfInterest": {"id": "urn:uuid:2341846f-af2a-487d-8ea0-e2161116dc35"}, "result": "urn:uuid:5c61153b-8293-b238-2543-cb6147b72f38", "procedure": {"id": "2024:aligner2024a", "implementedBy": "AlignerGeometryProcessor", "type": "sosa:Procedure", "ssn:hasInput": [{"id": "brdr:relevant_distance", "type": "ssn:Input", "input_value": {"type": "xsd:integer", "value": 1.4}}]}}, "observations": [{"type": "sosa:Observation", "has_feature_of_interest": "urn:uuid:863c0d38-4a8b-4184-a171-d6b0cd4ddc13", "made_by_sensor": "urn:uuid:501f1a88-5b1e-4dc9-aaa9-e42a260a3b33", "result_time": "2026-02-04T15:54:39", "id": "urn:uuid:8d078e26-eb19-4265-a393-515013600125", "observed_property": "brdr:area_overlap", "result": {"value": 764.67, "type": "float"}, "used_procedure": "brdr:observation_procedure_area_overlap", "used": "urn:uuid:5c61153b-8293-b238-2543-cb6147b72f38"}, {"type": "sosa:Observation", "has_feature_of_interest": "urn:uuid:863c0d38-4a8b-4184-a171-d6b0cd4ddc13", "made_by_sensor": "urn:uuid:501f1a88-5b1e-4dc9-aaa9-e42a260a3b33", "result_time": "2026-02-04T15:54:39", "id": "urn:uuid:382149fb-bcef-45b2-84d8-2f9df92d2570", "observed_property": "brdr:area_overlap_percentage", "result": {"value": 100, "type": "float"}, "used_procedure": "brdr:observation_procedure_area_overlap_percentage", "used": "urn:uuid:5c61153b-8293-b238-2543-cb6147b72f38"}, {"type": "sosa:Observation", "has_feature_of_interest": "urn:uuid:5c61153b-8293-b238-2543-cb6147b72f38", "made_by_sensor": "urn:uuid:501f1a88-5b1e-4dc9-aaa9-e42a260a3b33", "result_time": "2026-02-04T15:54:39", "id": "urn:uuid:218856a7-e912-47f5-aef8-f4169b363753", "observed_property": "brdr:area_overlap_full", "result": {"value": true, "type": "boolean"}, "used_procedure": "brdr:observation_procedure_area_overlap_full", "used": "urn:uuid:863c0d38-4a8b-4184-a171-d6b0cd4ddc13"}, {"type": "sosa:Observation", "has_feature_of_interest": "urn:uuid:5c61153b-8293-b238-2543-cb6147b72f38", "made_by_sensor": "urn:uuid:501f1a88-5b1e-4dc9-aaa9-e42a260a3b33", "result_time": "2026-02-04T15:54:39", "id": "urn:uuid:2be7c769-38dd-4fc7-b556-2e2332a088ec", "observed_property": "brdr:area", "result": {"value": 764.67, "type": "float"}, "used_procedure": "brdr:observation_procedure_area"}]}'

        self.metadata = json.loads(metadata_string)

    def test_empty_input(self):
        """Test if function returns empty dict when input is empty or observations are empty."""
        self.assertEqual(_reverse_metadata_observations_to_brdr_observation({}), {})
        self.assertEqual(
            _reverse_metadata_observations_to_brdr_observation({"obs": []}), {}
        )

    def test_happy_path_root_properties(self):
        """Test of observations (area, area_od) are placed correctly."""

        result = _reverse_metadata_observations_to_brdr_observation(self.metadata)

        self.assertEqual(result["area"], 764.67)

    def test_reference_features_aggregation(self):
        """Test if values of references are grouped correctly."""
        result = _reverse_metadata_observations_to_brdr_observation(self.metadata)
        capakey = "24126B0029/00S000"
        self.assertIn(capakey, result["reference_features"])
        self.assertEqual(result["reference_features"][capakey]["area"], 764.67)
        self.assertEqual(result["reference_features"][capakey]["percentage"], 100)

    def test_reverse_metadata_maps_reference_od(self):
        metadata = {
            "actuation": {"result": "urn:result", "reference_geometries": []},
            "observations": [
                {
                    "observed_property": "brdr:area_open_domain",
                    "result": {"value": 0.0, "type": "float"},
                    "has_feature_of_interest": "urn:result",
                }
            ],
        }
        result = _reverse_metadata_observations_to_brdr_observation(metadata)
        assert result["reference_od"] == {"area": 0.0}

    def test_reverse_metadata_maps_de9im_per_reference_feature(self):
        metadata = {
            "actuation": {
                "result": "urn:result",
                "reference_geometries": [
                    {"id": "urn:ref", "derived_from": {"id": "ref1"}}
                ],
            },
            "observations": [
                {
                    "observed_property": "brdr:de9im_intersection_matrix",
                    "result": {"value": "212101212", "type": "string"},
                    "has_feature_of_interest": "urn:ref",
                }
            ],
        }
        result = _reverse_metadata_observations_to_brdr_observation(metadata)
        assert result["reference_features"]["ref1"]["de9im"] == "212101212"


def test_metadata_observations_preserve_zero_and_false_values():
    process_result = {
        "observations": {
            "reference_features": {
                "ref1": {"area": 0.0, "percentage": 0.0, "de9im": "0FFFFFFF2"}
            },
            "full": False,
            "area": 0.0,
            "reference_od": {"area": 0.0},
        },
        "metadata": {"actuation": {"result": "urn:result"}},
    }
    obs = _get_metadata_observations_from_process_result(
        processResult=process_result, reference_lookup={"ref1": "urn:ref1"}
    )
    observed_props = {o["observed_property"]: o["result"]["value"] for o in obs}
    assert observed_props["brdr:de9im_intersection_matrix"] == "0FFFFFFF2"
    assert observed_props["brdr:area_overlap"] == 0.0
    assert observed_props["brdr:area_overlap_percentage"] == 0.0
    assert observed_props["brdr:area_overlap_full"] is False
    assert observed_props["brdr:area"] == 0.0
    assert observed_props["brdr:area_open_domain"] == 0.0


def test_compare_to_reference_point_uses_count_metric():
    aligner = Aligner()
    aligner.load_thematic_data(DictLoader({"theme_id_1": from_wkt("POINT (0 0)")}))
    aligner.load_reference_data(
        DictLoader({"ref_id_1": from_wkt("POLYGON ((-1 -1,-1 1,1 1,1 -1,-1 -1))")})
    )
    observation = aligner.compare_to_reference(from_wkt("POINT (0 0)"))
    assert observation["measure_type"] == "count"
    assert "ref_id_1" in observation["reference_features"]
    assert observation["reference_features"]["ref_id_1"]["count"] == 1.0
    assert observation["reference_features"]["ref_id_1"]["percentage"] == 100


def test_compare_to_reference_line_uses_length_metric():
    aligner = Aligner()
    aligner.load_thematic_data(
        DictLoader({"theme_id_1": from_wkt("LINESTRING (0 0,2 0)")})
    )
    aligner.load_reference_data(
        DictLoader({"ref_id_1": from_wkt("LINESTRING (1 0,3 0)")})
    )
    observation = aligner.compare_to_reference(from_wkt("LINESTRING (0 0,2 0)"))
    assert observation["measure_type"] == "length"
    assert "ref_id_1" in observation["reference_features"]
    assert observation["reference_features"]["ref_id_1"]["length"] == 1.0
    assert observation["reference_od"]["length"] == 1.0


def test_compare_to_reference_polygon_against_line_uses_line_metric():
    aligner = Aligner()
    aligner.load_thematic_data(
        DictLoader(
            {"theme_id_1": from_wkt("POLYGON ((0 0,0 2,2 2,2 0,0 0))")}
        )
    )
    aligner.load_reference_data(
        DictLoader({"ref_id_1": from_wkt("LINESTRING (-1 1, 3 1)")})
    )
    observation = aligner.compare_to_reference(
        from_wkt("POLYGON ((0 0,0 2,2 2,2 0,0 0))")
    )
    assert observation["reference_features"]["ref_id_1"]["measure_type"] == "length"
    assert observation["reference_features"]["ref_id_1"]["length"] == 2.0
    assert observation["reference_features"]["ref_id_1"]["percentage"] == 50.0


def test_compare_to_reference_line_against_polygon_uses_length_metric():
    aligner = Aligner()
    aligner.load_thematic_data(DictLoader({"theme_id_1": from_wkt("LINESTRING (0 0,4 0)")}))
    aligner.load_reference_data(
        DictLoader({"ref_id_1": from_wkt("POLYGON ((1 -1,1 1,3 1,3 -1,1 -1))")})
    )
    observation = aligner.compare_to_reference(from_wkt("LINESTRING (0 0,4 0)"))
    assert observation["reference_features"]["ref_id_1"]["measure_type"] == "length"
    assert observation["reference_features"]["ref_id_1"]["length"] == 2.0
    assert observation["reference_features"]["ref_id_1"]["percentage"] == 50.0


def test_compare_to_reference_adds_de9im_per_reference_feature():
    aligner = Aligner()
    aligner.load_thematic_data(
        DictLoader({"theme_id_1": from_wkt("POLYGON ((0 0,0 2,2 2,2 0,0 0))")})
    )
    aligner.load_reference_data(
        DictLoader({"ref_id_1": from_wkt("POLYGON ((1 0,1 2,3 2,3 0,1 0))")})
    )
    observation = aligner.compare_to_reference(
        from_wkt("POLYGON ((0 0,0 2,2 2,2 0,0 0))")
    )
    de9im = observation["reference_features"]["ref_id_1"]["de9im"]
    assert isinstance(de9im, str)
    assert len(de9im) == 9
