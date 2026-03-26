from brdr.aligner import Aligner
from brdr.configs import ProcessorConfig
from brdr.enums import SnapStrategy
from brdr.geometry_utils import geom_from_wkt
from brdr.loader import DictLoader
from brdr.processor import NetworkGeometryProcessor


def _run_network_with_strategy(
    strategy: SnapStrategy, thematic_wkt: str, relevant_distance: float
):
    processor_config = ProcessorConfig(
        snap_strategy=strategy,
        angle_threshold_degrees=150.0,
        snap_max_segment_length=2,
    )
    processor = NetworkGeometryProcessor(config=processor_config)
    aligner = Aligner(processor=processor, crs="EPSG:3812")

    thematic_id = "theme_line"
    thematic_data = {
        thematic_id: geom_from_wkt(thematic_wkt)
    }
    # Single reference line with open ends and a sharp bend.
    reference_data = {
        "ref_1": geom_from_wkt(
            "LINESTRING (0 0, 20 0, 40 0, 60 0, 70 20, 80 20)"
        ),
    }

    aligner.load_thematic_data(DictLoader(thematic_data))
    aligner.load_reference_data(DictLoader(reference_data))

    aligner_result = aligner.process(relevant_distances=[relevant_distance])
    process_results = aligner_result.get_results(aligner=aligner)
    return process_results[thematic_id][relevant_distance]["result"]


if __name__ == "__main__":
    """
    Compare NetworkGeometryProcessor output for:
    - SnapStrategy.PREFER_VERTICES
    - SnapStrategy.PREFER_VERTICES_ENDS_AND_ANGLES
    """
    # Simple single thematic line with endpoints near multiple candidate vertices.
    thematic_wkt = "LINESTRING (18 0, 60 20)"
    relevant_distance = 30.0

    result_prefer_vertices = _run_network_with_strategy(
        SnapStrategy.PREFER_VERTICES, thematic_wkt, relevant_distance
    )
    result_pref_vertices_ends_angles = _run_network_with_strategy(
        SnapStrategy.PREFER_VERTICES_ENDS_AND_ANGLES, thematic_wkt, relevant_distance
    )

    print("thematic_wkt:", thematic_wkt)
    print("relevant_distance:", relevant_distance)
    print()
    print("PREFER_VERTICES result:")
    print(result_prefer_vertices.wkt)
    print()
    print("PREFER_VERTICES_ENDS_AND_ANGLES result:")
    print(result_pref_vertices_ends_angles.wkt)
    print()
    print("Equal:", result_prefer_vertices.equals(result_pref_vertices_ends_angles))
