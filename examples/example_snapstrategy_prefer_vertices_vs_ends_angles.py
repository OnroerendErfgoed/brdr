from shapely.geometry import LineString, MultiLineString, Point

from brdr.enums import SnapStrategy
from brdr.geometry_utils import snap_geometry_to_reference


def _snap_single_point(point: Point, reference, strategy: SnapStrategy, tolerance: float):
    snapped = snap_geometry_to_reference(
        point,
        reference,
        snap_strategy=strategy,
        tolerance=tolerance,
        angle_threshold_degrees=150.0,
    )
    return tuple(snapped.coords[0])


if __name__ == "__main__":
    """
    Demonstrates the difference between:
    - SnapStrategy.PREFER_VERTICES
    - SnapStrategy.PREFER_VERTICES_ENDS_AND_ANGLES
    """
    reference = MultiLineString(
        [
            LineString(
                [
                    (0, 0),   # end vertex
                    (2, 0),   # regular vertex
                    (4, 0),   # regular vertex
                    (6, 0),   # angle vertex (90 deg turn)
                    (6, 2),   # angle vertex
                    (8, 2),   # end vertex
                ]
            )
        ]
    )

    tolerance = 2.5

    # Case 1: end vertex is farther than a regular vertex but still within tolerance.
    # PREFER_VERTICES picks closest vertex (2,0).
    # PREFER_VERTICES_ENDS_AND_ANGLES prefers end vertex first, so picks (0,0).
    case_end_priority = Point(1.8, 0.1)

    # Case 2: no end vertex nearby; angle vertex and regular vertex are both within tolerance.
    # PREFER_VERTICES picks closest regular vertex (4,0).
    # PREFER_VERTICES_ENDS_AND_ANGLES prefers angle vertex over regular, so picks (6,0).
    case_angle_priority = Point(4.1, 0.2)

    test_points = [
        ("case_end_priority", case_end_priority),
        ("case_angle_priority", case_angle_priority),
    ]

    print("Reference vertices:", list(reference.geoms[0].coords))
    print("Tolerance:", tolerance)
    print()

    for label, pt in test_points:
        snapped_prefer_vertices = _snap_single_point(
            pt, reference, SnapStrategy.PREFER_VERTICES, tolerance
        )
        snapped_ends_angles = _snap_single_point(
            pt,
            reference,
            SnapStrategy.PREFER_VERTICES_ENDS_AND_ANGLES,
            tolerance,
        )

        print(label)
        print("  input:", tuple(pt.coords[0]))
        print("  PREFER_VERTICES:", snapped_prefer_vertices)
        print("  PREFER_VERTICES_ENDS_AND_ANGLES:", snapped_ends_angles)
        print()
