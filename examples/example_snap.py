import math

from shapely.geometry import Point, LineString


def distance(p1, p2):
    return math.hypot(p1[0] - p2[0], p1[1] - p2[1])


def snap_point_to_nodes(point, target_lines, tolerance):
    closest_point = None
    min_dist = float("inf")
    for line in target_lines:
        for coord in line.coords:
            d = distance(point, coord)
            if d < min_dist and d <= tolerance:
                min_dist = d
                closest_point = coord
    return closest_point


def snap_point_to_edges(point, target_lines, tolerance):
    p = Point(point)
    closest_proj = None
    min_dist = float("inf")
    for line in target_lines:
        proj = line.interpolate(line.project(p))
        d = p.distance(proj)
        if d < min_dist and d <= tolerance:
            min_dist = d
            closest_proj = (proj.x, proj.y)
    return closest_proj


def is_parallel(line1, line2, angle_tolerance=5):
    def angle(line):
        x1, y1 = line.coords[0]
        x2, y2 = line.coords[-1]
        return math.degrees(math.atan2(y2 - y1, x2 - x1))

    a1 = angle(line1)
    a2 = angle(line2)
    diff = abs((a1 - a2 + 180) % 360 - 180)
    return diff <= angle_tolerance


def snap_geometry(source_line, target_lines, tolerance):
    snapped_coords = []
    for coord in source_line.coords:
        # Try node snap
        snapped = snap_point_to_nodes(coord, target_lines, tolerance)
        if snapped:
            snapped_coords.append(snapped)
            continue

        # Try edge snap
        snapped = snap_point_to_edges(coord, target_lines, tolerance)
        if snapped:
            snapped_coords.append(snapped)
            continue

        # No snap
        snapped_coords.append(coord)

    return LineString(snapped_coords)


# Example usage
if __name__ == "__main__":
    source = LineString([(0, 0), (5, 5), (10, 10)])
    target1 = LineString([(0, 1), (10, 11)])
    target2 = LineString([(5, 0), (5, 10)])
    tolerance = 1.5

    snapped = snap_geometry(source, [target1, target2], tolerance)

    print("Origineel:", list(source.coords))
    print("Gesnapt:  ", list(snapped.coords))
