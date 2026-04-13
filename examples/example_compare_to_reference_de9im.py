from brdr.aligner import Aligner
from brdr.loader import DictLoader
from shapely import from_wkt


if __name__ == "__main__":
    """Example: inspect generic observation metrics and DE-9IM per reference feature."""

    thematic = {"line_1": from_wkt("LINESTRING (0 0, 10 0)")}
    reference = {
        "poly_ref": from_wkt("POLYGON ((2 -1, 8 -1, 8 1, 2 1, 2 -1))"),
        "line_ref": from_wkt("LINESTRING (0 0, 10 0)"),
        "point_ref": from_wkt("POINT (5 0)"),
    }

    aligner = Aligner()
    aligner.load_thematic_data(DictLoader(thematic))
    aligner.load_reference_data(DictLoader(reference))

    observation = aligner.compare_to_reference(thematic["line_1"])

    print("measure_type:", observation["measure_type"])
    print("reference_od:", observation["reference_od"])
    print("reference features:")
    for ref_id, info in observation["reference_features"].items():
        print(
            f"- {ref_id}: "
            f"measure_type={info.get('measure_type')}, "
            f"percentage={info.get('percentage')}, "
            f"de9im={info.get('de9im')}"
        )
