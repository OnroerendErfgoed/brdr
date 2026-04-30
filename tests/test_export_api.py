import json
from types import SimpleNamespace

from shapely import from_wkt

from brdr.aligner import AlignerResult
from brdr.loader import DictLoader


def _build_feature_collection():
    source_dict = {
        "a": from_wkt("LINESTRING (0 0, 1 1)"),
        "b": from_wkt("POINT (2 3)"),
    }
    source_props = {
        "a": {"name": "alpha", "weight": 1},
        "b": {"name": "beta", "weight": 2},
    }
    return DictLoader(
        data_dict=source_dict,
        data_dict_properties=source_props,
        is_reference=False,
    ).load_data()


def test_feature_collection_export_profiles_and_shortcuts(tmp_path):
    collection = _build_feature_collection()

    gdf_full = collection.to_gdf(profile="full")
    assert collection.id_fieldname in gdf_full.columns
    assert "geometry" in gdf_full.columns
    assert "name" in gdf_full.columns

    gdf_minimal = collection.to_gdf(profile="minimal")
    assert set(gdf_minimal.columns) == {collection.id_fieldname, "geometry"}

    geojson_payload = collection.to_geojson_dict(profile="minimal")
    assert geojson_payload["type"] == "FeatureCollection"

    export_gdf = collection.export("gdf", profile="minimal")
    assert set(export_gdf.columns) == {collection.id_fieldname, "geometry"}

    out_geojson = tmp_path / "collection.geojson"
    info = collection.export("geojson", path=str(out_geojson), profile="full")
    assert info["format"] == "geojson"
    assert out_geojson.exists()
    with out_geojson.open("r", encoding="utf-8") as fp:
        payload = json.load(fp)
    assert payload["type"] == "FeatureCollection"


def test_aligner_result_export_shortcuts_smoke():
    thematic = _build_feature_collection()
    aligner_stub = SimpleNamespace(thematic_data=thematic, crs=thematic.crs)
    process_results = {
        "a": {
            1.0: {
                "result": thematic.features["a"].geometry,
                "result_diff": thematic.features["a"].geometry,
                "properties": {"score": 1.0},
            }
        }
    }
    result = AlignerResult(process_results=process_results)

    gdf = result.to_gdf(aligner_stub, profile="minimal")
    assert thematic.id_fieldname in gdf.columns
    assert "result" in gdf.columns or "geometry" in gdf.columns

    geojson_payload = result.to_geojson(aligner_stub, profile="minimal")
    assert geojson_payload["type"] == "FeatureCollection"

    exported = result.export(aligner_stub, format="gdf", profile="full")
    assert thematic.id_fieldname in exported.columns
