from shapely import from_wkt

from brdr.aligner import Aligner
from brdr.descriptor import AlignerDescriptor
from brdr.loader import DictLoader


def test_geometry_descriptor_get_actual_observation_for_line():
    aligner = Aligner()
    aligner.load_thematic_data(
        DictLoader({"theme_1": from_wkt("LINESTRING (0 0,2 0)")})
    )
    aligner.load_reference_data(
        DictLoader({"ref_1": from_wkt("LINESTRING (1 0,3 0)")})
    )
    descriptor = AlignerDescriptor()
    process_result = {"result": from_wkt("LINESTRING (0 0,2 0)"), "properties": {}}

    obs = descriptor.get_actual_observation(
        aligner=aligner, process_result=process_result, cache_key="line1"
    )
    assert obs is not None
    assert obs["measure_type"] == "length"
    assert "reference_features" in obs


def test_geometry_descriptor_get_base_observation_from_brdr_payload():
    descriptor = AlignerDescriptor()
    payload = {
        "alignment_date": "2026-01-01T00:00:00",
        "brdr_version": "x",
        "reference_source": "s",
        "full": True,
        "area": 1.0,
        "length": 0.0,
        "count": 0.0,
        "measure_type": "area",
        "reference_features": {},
        "reference_od": None,
    }
    obs = descriptor.get_base_observation(
        feature_properties={"meta": payload}, metadata_field="meta", cache_key="base1"
    )
    assert obs is not None
    assert obs["measure_type"] == "area"
    assert obs["area"] == 1.0
