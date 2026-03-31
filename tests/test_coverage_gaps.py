import re

import pytest
from shapely import from_wkt

from brdr.aligner import Aligner
from brdr.configs import ProcessorConfig
from brdr.constants import EVALUATION_FIELD_NAME
from brdr.constants import VERSION_DATE
from brdr.enums import Evaluation
from brdr.enums import OpenDomainStrategy
from brdr.evaluator import AlignerEvaluator
from brdr.evaluator import _reverse_metadata_observations_to_brdr_observation
from brdr.loader import DictLoader
from brdr.loader import GeoDataFrameLoader
from brdr.loader import OGCFeatureAPIReferenceLoader
from brdr.loader import WFSReferenceLoader
from brdr.metadata import get_metadata_observations_from_process_result
from brdr.metadata import reverse_metadata_observations_to_brdr_observation
from brdr.processor import SnapGeometryProcessor


def _sample_feature_collection(fid_name: str = "OIDN", fid_value: str = "ref1"):
    return {
        "type": "FeatureCollection",
        "features": [
            {
                "type": "Feature",
                "properties": {fid_name: fid_value},
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [[[0, 0], [0, 2], [2, 2], [2, 0], [0, 0]]],
                },
            }
        ],
    }


def test_evaluate_marks_non_selected_thematic_ids_as_not_evaluated():
    aligner = Aligner()
    aligner.load_thematic_data(
        DictLoader(
            {
                "theme_1": from_wkt("POLYGON ((0 0,0 2,2 2,2 0,0 0))"),
                "theme_2": from_wkt("POLYGON ((3 0,3 2,5 2,5 0,3 0))"),
            }
        )
    )
    aligner.load_reference_data(
        DictLoader({"ref_1": from_wkt("POLYGON ((0 0,0 2,2 2,2 0,0 0))")})
    )

    result = aligner.evaluate(relevant_distances=[0], thematic_ids=["theme_1"])
    eval_theme_2 = result.results["theme_2"][0]["properties"][EVALUATION_FIELD_NAME]
    assert eval_theme_2 == Evaluation.NOT_EVALUATED


def test_reverse_metadata_observation_wrapper_returns_empty_for_invalid_payload():
    assert _reverse_metadata_observations_to_brdr_observation({}) == {}


def test_snap_processor_od_polygon_reference_changes_result_by_strategy():
    thematic = {"theme_id": from_wkt("POLYGON ((0 0,0 2,10 2,10 0,0 0))")}
    reference = {
        "ref_id": from_wkt("POLYGON ((2 -1,2 3,8 3,8 -1,2 -1))"),
    }

    aligner_exclude = Aligner(
        processor=SnapGeometryProcessor(
            config=ProcessorConfig(od_strategy=OpenDomainStrategy.EXCLUDE)
        )
    )
    aligner_exclude.load_thematic_data(DictLoader(thematic))
    aligner_exclude.load_reference_data(DictLoader(reference))
    result_exclude = aligner_exclude.process([0.5]).results["theme_id"][0.5]["result"]

    aligner_as_is = Aligner(
        processor=SnapGeometryProcessor(
            config=ProcessorConfig(od_strategy=OpenDomainStrategy.AS_IS)
        )
    )
    aligner_as_is.load_thematic_data(DictLoader(thematic))
    aligner_as_is.load_reference_data(DictLoader(reference))
    result_as_is = aligner_as_is.process([0.5]).results["theme_id"][0.5]["result"]

    assert result_as_is.area > result_exclude.area


def test_snap_processor_od_ignores_non_polygon_reference():
    thematic = {"theme_id": from_wkt("LINESTRING (0 0,10 0)")}
    reference = {"ref_id": from_wkt("LINESTRING (2 0,8 0)")}

    aligner_exclude = Aligner(
        processor=SnapGeometryProcessor(
            config=ProcessorConfig(od_strategy=OpenDomainStrategy.EXCLUDE)
        )
    )
    aligner_exclude.load_thematic_data(DictLoader(thematic))
    aligner_exclude.load_reference_data(DictLoader(reference))
    result_exclude = aligner_exclude.process([0.5]).results["theme_id"][0.5]["result"]

    aligner_as_is = Aligner(
        processor=SnapGeometryProcessor(
            config=ProcessorConfig(od_strategy=OpenDomainStrategy.AS_IS)
        )
    )
    aligner_as_is.load_thematic_data(DictLoader(thematic))
    aligner_as_is.load_reference_data(DictLoader(reference))
    result_as_is = aligner_as_is.process([0.5]).results["theme_id"][0.5]["result"]

    assert result_as_is.length == result_exclude.length


def test_loader_swallow_invalid_version_date():
    loader = DictLoader(
        {"id1": from_wkt("POLYGON ((0 0,0 1,1 1,1 0,0 0))")},
        data_dict_properties={"id1": {"d": "not-a-date"}},
    )
    loader.versiondate_info = {"name": "d", "format": "%Y-%m-%d"}
    data = loader.load_data()
    assert VERSION_DATE not in data.features["id1"].properties


def test_geodataframe_loader_missing_id_property_raises():
    gpd = pytest.importorskip("geopandas")
    gdf = gpd.GeoDataFrame.from_features(_sample_feature_collection()["features"])
    loader = GeoDataFrameLoader(_input=gdf, id_property="missing_id")
    with pytest.raises(KeyError):
        loader.load_data()


def test_ogc_feature_api_loader_happy_path(monkeypatch, requests_mock):
    aligner = Aligner()
    aligner.load_thematic_data(
        DictLoader({"theme_id": from_wkt("POLYGON ((0 0,0 2,2 2,2 0,0 0))")})
    )

    base_url = "https://example.com/ogc"
    requests_mock.add(
        requests_mock.GET,
        re.compile(r"^https://example\.com/ogc/collections$"),
        json={"collections": [{"id": "ADP"}]},
        status=200,
    )
    requests_mock.add(
        requests_mock.GET,
        re.compile(r"^https://example\.com/ogc/collections/ADP$"),
        json={"crs": ["http://www.opengis.net/def/crs/EPSG/0/31370"]},
        status=200,
    )
    monkeypatch.setattr(
        "brdr.loader.get_collection_by_partition",
        lambda **kwargs: _sample_feature_collection(),
    )

    loader = OGCFeatureAPIReferenceLoader(
        url=base_url,
        id_property="OIDN",
        collection="ADP",
        aligner=aligner,
        partition=500,
    )
    loaded = loader.load_data()
    assert "ref1" in loaded.features
    assert VERSION_DATE in loaded.source


def test_wfs_loader_happy_path(monkeypatch, requests_mock):
    aligner = Aligner()
    aligner.load_thematic_data(
        DictLoader({"theme_id": from_wkt("POLYGON ((0 0,0 2,2 2,2 0,0 0))")})
    )
    capabilities_xml = """<?xml version="1.0" encoding="UTF-8"?>
    <WFS_Capabilities xmlns:wfs="http://www.opengis.net/wfs/2.0">
      <FeatureTypeList>
        <FeatureType>
          <Name>ns:layer</Name>
          <DefaultCRS>urn:ogc:def:crs:EPSG::31370</DefaultCRS>
        </FeatureType>
      </FeatureTypeList>
    </WFS_Capabilities>
    """
    requests_mock.add(
        requests_mock.GET,
        re.compile(r"^https://example\.com/wfs.*"),
        body=capabilities_xml,
        status=200,
        content_type="application/xml",
    )
    monkeypatch.setattr(
        "brdr.loader.get_collection_by_partition",
        lambda **kwargs: _sample_feature_collection(),
    )

    loader = WFSReferenceLoader(
        url="https://example.com/wfs",
        id_property="OIDN",
        typename="ns:layer",
        aligner=aligner,
    )
    loaded = loader.load_data()
    assert "ref1" in loaded.features


def test_wfs_loader_missing_typename_raises(requests_mock):
    aligner = Aligner()
    aligner.load_thematic_data(
        DictLoader({"theme_id": from_wkt("POLYGON ((0 0,0 2,2 2,2 0,0 0))")})
    )
    capabilities_xml = """<?xml version="1.0" encoding="UTF-8"?>
    <WFS_Capabilities xmlns:wfs="http://www.opengis.net/wfs/2.0">
      <FeatureTypeList>
        <FeatureType>
          <Name>ns:other</Name>
          <DefaultCRS>urn:ogc:def:crs:EPSG::31370</DefaultCRS>
        </FeatureType>
      </FeatureTypeList>
    </WFS_Capabilities>
    """
    requests_mock.add(
        requests_mock.GET,
        re.compile(r"^https://example\.com/wfs.*"),
        body=capabilities_xml,
        status=200,
        content_type="application/xml",
    )

    loader = WFSReferenceLoader(
        url="https://example.com/wfs",
        id_property="OIDN",
        typename="ns:layer",
        aligner=aligner,
    )
    with pytest.raises(ValueError, match="Typename"):
        loader.load_data()


def test_metadata_observations_for_length_metric_and_area_od_fallback():
    process_result = {
        "observations": {
            "measure_type": "length",
            "reference_features": {
                "ref1": {
                    "measure_type": "length",
                    "length": 4.0,
                    "percentage": 50.0,
                }
            },
            "full": False,
            "length": 8.0,
            "area_od": {"length": 1.5},
        },
        "metadata": {"actuation": {"result": "urn:result"}},
    }
    observations = get_metadata_observations_from_process_result(
        processResult=process_result,
        reference_lookup={"ref1": "urn:ref1"},
    )
    by_prop = {o["observed_property"]: o["result"]["value"] for o in observations}
    assert by_prop["brdr:length_overlap"] == 4.0
    assert by_prop["brdr:length_overlap_percentage"] == 50.0
    assert by_prop["brdr:length"] == 8.0
    assert by_prop["brdr:length_open_domain"] == 1.5


def test_reverse_metadata_supports_length_and_count_observations():
    metadata = {
        "actuation": {
            "result": "urn:result",
            "reference_geometries": [{"id": "urn:ref", "derived_from": {"id": "r1"}}],
        },
        "observations": [
            {
                "observed_property": "brdr:length",
                "result": {"value": 10.0, "type": "float"},
                "has_feature_of_interest": "urn:result",
            },
            {
                "observed_property": "brdr:length_open_domain",
                "result": {"value": 2.0, "type": "float"},
                "has_feature_of_interest": "urn:result",
            },
            {
                "observed_property": "brdr:length_overlap",
                "result": {"value": 8.0, "type": "float"},
                "has_feature_of_interest": "urn:ref",
            },
            {
                "observed_property": "brdr:length_overlap_percentage",
                "result": {"value": 100.0, "type": "float"},
                "has_feature_of_interest": "urn:ref",
            },
        ],
    }
    result = reverse_metadata_observations_to_brdr_observation(metadata)
    assert result["measure_type"] == "length"
    assert result["reference_od"] == {"length": 2.0}
    assert result["reference_features"]["r1"]["length"] == 8.0
    assert result["reference_features"]["r1"]["full"] is True


def _valid_brdr_obs(
    *,
    full: bool,
    reference_features: dict,
    reference_od=None,
    measure_type: str = "area",
):
    return {
        "alignment_date": None,
        "brdr_version": None,
        "reference_source": None,
        "full": full,
        "area": 0.0,
        "length": None,
        "count": None,
        "measure_type": measure_type,
        "reference_features": reference_features,
        "reference_od": reference_od,
    }


def test_evaluator_comparison_properties_support_metric_fallbacks():
    evaluator = AlignerEvaluator()
    aligner = Aligner()
    process_result = {
        "result": from_wkt("LINESTRING (0 0, 2 0)"),
        "observations": _valid_brdr_obs(
            full=False,
            measure_type="invalid",
            reference_features={"r1": {"full": True, "area": 2.0}},
            reference_od=None,
        ),
    }
    base_obs = _valid_brdr_obs(
        full=False,
        measure_type="count",
        reference_features={"r1": {"full": True, "count": 2.0, "measure_type": "count"}},
        reference_od=None,
    )

    props = evaluator.get_observation_comparison_properties(
        aligner=aligner,
        process_result=process_result,
        base_brdr_observation=base_obs,
    )
    assert props[EVALUATION_FIELD_NAME] == Evaluation.EQUALITY_BY_ID
    assert props["brdr_diff_area"] == 0.0
    assert props["brdr_diff_percentage"] == 0.0


def test_evaluator_comparison_properties_with_zero_od_marks_full_equality():
    evaluator = AlignerEvaluator()
    aligner = Aligner()
    process_result = {
        "result": from_wkt("POLYGON ((0 0,0 2,2 2,2 0,0 0))"),
        "observations": _valid_brdr_obs(
            full=True,
            reference_features={"r1": {"full": True, "area": 4.0}},
            reference_od={"area": 0.0},
        ),
    }
    base_obs = _valid_brdr_obs(
        full=True,
        reference_features={"r1": {"full": True, "area": 4.0}},
        reference_od={"area": 0.0},
    )

    props = evaluator.get_observation_comparison_properties(
        aligner=aligner,
        process_result=process_result,
        base_brdr_observation=base_obs,
    )
    assert props[EVALUATION_FIELD_NAME] == Evaluation.EQUALITY_BY_ID_AND_FULL_REFERENCE
    assert props["brdr_od_alike"] is True
