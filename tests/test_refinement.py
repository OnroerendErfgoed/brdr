from shapely.geometry import LineString, Point, MultiPoint
from shapely import STRtree

from brdr.geometry_utils import (
    _refine_line_near_reference_vertices,
    _refine_multipolygon_near_reference_vertices,
    geom_from_wkt,
    snap_geometry_to_reference,
)


def test_refine_line_near_reference_vertices_is_local():
    line = LineString([(0.0, 0.0), (100.0, 0.0)])
    ref_points = [Point(50.0, 0.2)]
    tree = STRtree(ref_points)

    refined = _refine_line_near_reference_vertices(
        line,
        ref_vertex_points=ref_points,
        ref_vertex_tree=tree,
        tolerance=1.0,
        max_segment_length=2.0,
    )

    # Local refinement should add only a few points around projection (~x=50)
    assert len(refined.coords) <= 6
    xs = [c[0] for c in refined.coords]
    assert any(abs(x - 50.0) < 1e-6 for x in xs)


def test_snap_multipoint_ignores_max_segment_length():
    points = MultiPoint([(0.0, 0.0), (1.0, 1.0)])
    ref = MultiPoint([(0.1, 0.0), (1.0, 1.2)])
    snapped = snap_geometry_to_reference(
        points,
        reference=ref,
        max_segment_length=0.1,
        tolerance=2.0,
    )
    assert snapped.geom_type == "MultiPoint"
    assert len(list(snapped.geoms)) == 2


def test_refine_multipolygon_near_reference_vertices_is_local():
    poly = geom_from_wkt("MULTIPOLYGON (((0 0, 100 0, 100 10, 0 10, 0 0)))")
    ref = geom_from_wkt("MULTIPOINT ((50 0.2), (50 9.8))")
    refined = _refine_multipolygon_near_reference_vertices(
        poly,
        reference=ref,
        tolerance=1.0,
        max_segment_length=2.0,
    )
    ring = list(list(refined.geoms)[0].exterior.coords)
    # Should be refined locally around x~50; not globally densified.
    assert len(ring) < 20
