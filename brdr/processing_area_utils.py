from typing import Callable

from shapely import GeometryCollection, MultiPolygon, Point, Polygon, line_merge, make_valid
from shapely.geometry import LineString, LinearRing
from shapely.geometry.base import BaseGeometry
from shapely.ops import nearest_points, snap

from brdr.geometry_utils import safe_difference, safe_intersection, safe_unary_union


def align_polygon_boundary_in_processing_area(
    *,
    geometry: BaseGeometry,
    processing_area: BaseGeometry,
    relevant_distance: float,
    correction_distance: float,
    align_lines_fn: Callable[[BaseGeometry], BaseGeometry | None],
) -> BaseGeometry | None:
    """
    Replace polygon boundary parts inside `processing_area` by aligned line parts.

    Steps:
    1. Cut boundary lines by `processing_area`.
    2. Keep cut endpoints (anchors).
    3. Align each in-scope line part with `align_lines_fn` (Network processor callback).
    4. Reconnect original anchors to the aligned part endpoints/nearest points.
    5. Merge with out-of-scope boundary and rebuild polygon rings.
    """

    if geometry is None or geometry.is_empty:
        return geometry
    if not isinstance(geometry, (Polygon, MultiPolygon)):
        return None
    if processing_area is None or processing_area.is_empty:
        return geometry

    processing_area_valid = make_valid(processing_area)
    if processing_area_valid is None or processing_area_valid.is_empty:
        return geometry

    boundary = geometry.boundary
    boundary_in_scope = make_valid(safe_intersection(boundary, processing_area_valid))
    if boundary_in_scope is None or boundary_in_scope.is_empty:
        return geometry
    segment_parts = [
        g
        for g in getattr(boundary_in_scope, "geoms", [boundary_in_scope])
        if isinstance(g, LineString) and not g.is_empty and len(g.coords) >= 2
    ]
    if not segment_parts:
        return geometry

    aligned_segments = []
    connector_lines = []
    endpoint_snap_tol = max(
        float(correction_distance) * 5.0, float(relevant_distance) * 0.05
    )
    anchor_connect_max = max(float(correction_distance) * 50.0, float(relevant_distance) * 3.0)

    for seg in segment_parts:
        # Keep intersection anchors from the original cut line.
        start_anchor = Point(seg.coords[0])
        end_anchor = Point(seg.coords[-1])

        # Align in-scope line part with network-based callback.
        aligned_seg = align_lines_fn(seg)
        if aligned_seg is None or aligned_seg.is_empty:
            aligned_seg = seg
        aligned_seg = make_valid(aligned_seg)
        if aligned_seg is None or aligned_seg.is_empty:
            aligned_seg = seg

        # Keep only the aligned section inside the processing area.
        aligned_in_scope = make_valid(safe_intersection(aligned_seg, processing_area_valid))
        if aligned_in_scope is None or aligned_in_scope.is_empty:
            aligned_in_scope = seg

        # Snap toward anchors to favor exact continuity at cut points.
        anchor_points = [start_anchor, end_anchor]
        for anchor in anchor_points:
            aligned_in_scope = snap(aligned_in_scope, anchor, endpoint_snap_tol)
        aligned_in_scope = make_valid(aligned_in_scope)
        if aligned_in_scope is None or aligned_in_scope.is_empty:
            aligned_in_scope = seg

        # Explicitly connect original cut anchors to aligned result.
        for anchor in anchor_points:
            try:
                _, nearest_aligned = nearest_points(anchor, aligned_in_scope)
            except Exception:
                continue
            if nearest_aligned is None or nearest_aligned.is_empty:
                continue
            dist = anchor.distance(nearest_aligned)
            if dist <= endpoint_snap_tol:
                continue
            if dist > anchor_connect_max:
                continue
            connector_lines.append(
                LineString(
                    [
                        (float(anchor.x), float(anchor.y)),
                        (float(nearest_aligned.x), float(nearest_aligned.y)),
                    ]
                )
            )
        aligned_segments.append(aligned_in_scope)

    if not aligned_segments:
        return geometry

    aligned_in_scope_all = make_valid(safe_unary_union(aligned_segments))
    if connector_lines:
        aligned_in_scope_all = make_valid(
            safe_unary_union([aligned_in_scope_all, safe_unary_union(connector_lines)])
        )
    if aligned_in_scope_all is None or aligned_in_scope_all.is_empty:
        return geometry

    boundary_out_scope = make_valid(safe_difference(boundary, processing_area_valid))
    merged_boundary = make_valid(
        safe_unary_union([boundary_out_scope, aligned_in_scope_all])
    )
    if merged_boundary is None or merged_boundary.is_empty:
        return geometry

    # Stabilize merged boundary against the original linework.
    merge_snap_tol = max(float(correction_distance) * 10.0, float(relevant_distance) * 0.05)
    merged_boundary = snap(merged_boundary, boundary, merge_snap_tol)
    merged_boundary = make_valid(merged_boundary)
    if merged_boundary is None or merged_boundary.is_empty:
        return geometry

    merged_lines = line_merge(merged_boundary)
    ring_candidates = []
    for g in getattr(merged_lines, "geoms", [merged_lines]):
        if not isinstance(g, LineString):
            continue
        if g.is_empty or len(g.coords) < 4:
            continue
        if not g.is_ring:
            continue
        try:
            ring = LinearRing(g.coords)
            poly = Polygon(ring)
        except Exception:
            continue
        if poly.is_empty:
            continue
        ring_candidates.append(poly)
    if not ring_candidates:
        return geometry

    kept = []
    for poly in ring_candidates:
        overlap_area = poly.intersection(geometry).area
        if overlap_area > 0:
            kept.append(poly)
    if not kept:
        return geometry

    rebuilt = make_valid(safe_unary_union(kept))
    if rebuilt is None or rebuilt.is_empty:
        return geometry

    # Strictly apply edits inside processing area only.
    outside_original = safe_difference(geometry, processing_area_valid)
    inside_rebuilt = safe_intersection(rebuilt, processing_area_valid)
    merged_result = make_valid(safe_unary_union([outside_original, inside_rebuilt]))
    if merged_result is None:
        return geometry
    return merged_result
