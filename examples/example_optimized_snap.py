from shapely.geometry import LineString, Polygon, MultiLineString, MultiPoint, Point
from shapely.lib import segmentize
from shapely.ops import nearest_points, unary_union
from shapely.validation import make_valid

from brdr.enums import SnapStrategy
from brdr.geometry_utils import get_coords_from_geometry, closest_line, _get_line_substring, to_multi, \
    _snapped_point_by_snapstrategy, buffer_neg_pos


def snap_geometry_to_reference(geometries, reference, snap_strategy, tolerance):
    reference_coords = list(get_coords_from_geometry(reference))
    ref_lines = get_reference_lines(reference)
    lines = []

    for geom in geometries:
        coords = list(geom.coords)
        coordinates = []

        for idx, coord in enumerate(coords):
            if idx == 0:
                continue
            p_start = Point(coords[idx - 1])
            p_end = Point(coords[idx])

            p_start_snapped, bool_start_snapped = get_snapped_point(p_start, ref_lines, reference_coords, snap_strategy, tolerance)
            p_end_snapped, bool_end_snapped = get_snapped_point(p_end, ref_lines, reference_coords, snap_strategy, tolerance)

            coordinates.append(p_start_snapped.coords[0])

            if not bool_start_snapped or not bool_end_snapped:
                coordinates.append(p_end_snapped.coords[0])
                continue

            reference_line, distance = closest_line(ref_lines, p_end)
            distance_start_end = p_start.distance(p_end)
            line_substring = _get_line_substring(reference_line, p_start_snapped, p_end_snapped, distance_start_end)

            for p in line_substring.coords:
                point = Point(p)
                if (point.distance(p_start) + point.distance(p_end)) / 2 <= tolerance:
                    coordinates.append(p)
            coordinates.append(p_end_snapped.coords[0])

        lines.append(make_valid(LineString(coordinates)))

    return unary_union(lines)

def get_reference_lines(reference):
    reference = to_multi(reference, geomtype="Polygon")
    ref_lines = []
    for g in reference.geoms:
        ref_lines.append(g.exterior)
        ref_lines.extend([interior for interior in g.interiors])
    return unary_union(ref_lines)

def get_snapped_point(point, ref_lines, reference_coords, snap_strategy, tolerance):
    if not ref_lines.is_empty:
        p1, p2 = nearest_points(point, ref_lines)
    else:
        p1, p2 = None, None

    if len(reference_coords) != 0:
        p1_vertices, p2_vertices = nearest_points(point, MultiPoint(reference_coords))
    else:
        p1_vertices, p2_vertices = None, None

    return _snapped_point_by_snapstrategy(point, p2, p2_vertices, snap_strategy, tolerance)

def snap_line_to_polygon(geometry, reference, snap_strategy=SnapStrategy.NO_PREFERENCE, max_segment_length=-1, tolerance=1):
    if geometry is None or geometry.is_empty or reference is None or reference.is_empty:
        return geometry
    if max_segment_length > 0:
        geometry = segmentize(geometry, max_segment_length=max_segment_length)
    if geometry.geom_type == "LineString":
        geometry = MultiLineString([geometry])
    geoms =[geom for geom in geometry.geoms]
    return snap_geometry_to_reference(geoms, reference, snap_strategy, tolerance)

def snap_polygon_to_polygon(geometry, reference, snap_strategy=SnapStrategy.PREFER_VERTICES, max_segment_length=-1, tolerance=1, correction_distance=0.01):
    if geometry is None or geometry.is_empty or reference is None or reference.is_empty or not reference.geom_type in ("Polygon", "MultiPolygon"):
        return geometry
    if max_segment_length > 0:
        geometry = segmentize(geometry, max_segment_length=max_segment_length)
    geometry = to_multi(geometry, geomtype="Polygon")
    geoms = [geom.exterior for geom in geometry.geoms]
    result = snap_geometry_to_reference(geoms, reference, snap_strategy, tolerance)
    return buffer_neg_pos(result, correction_distance)

# Voorbeeld gebruik
line = LineString([(0, 0), (1, 1), (2, 2)])
polygon1 = Polygon([(0, 0), (0, 2), (2, 2), (2, 0)])
polygon2 = Polygon([(1, 1), (1, 3), (3, 3), (3, 1)])

tolerance = 0.5

snapped_line = snap_line_to_polygon(line, polygon1, tolerance=tolerance)
snapped_polygon = snap_polygon_to_polygon(polygon1, polygon2, tolerance=tolerance)

print(f"Gesnapte lijn: {snapped_line}")
print(f"Gesnapte polygon: {snapped_polygon}")