import json
import unittest

from shapely import from_wkt

from brdr.aligner import Aligner, _reverse_metadata_observations_to_brdr_observation
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
