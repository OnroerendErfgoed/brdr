import logging
from collections import Counter
from itertools import combinations, islice
from math import pi, inf, dist

import numpy as np
from shapely import (
    GEOSException,
    equals,
    shortest_line,
    GeometryCollection,
    MultiLineString,
    MultiPolygon,
    Polygon,
    MultiPoint,
)
from shapely import STRtree
from shapely import buffer
from shapely import difference
from shapely import from_wkt
from shapely import get_exterior_ring
from shapely import get_interior_ring
from shapely import get_num_interior_rings
from shapely import get_parts
from shapely import intersection
from shapely import is_empty
from shapely import make_valid
from shapely import polygons
from shapely import segmentize
from shapely import snap
from shapely import symmetric_difference
from shapely import to_wkt
from shapely import unary_union
from shapely import union
from shapely.geometry.base import BaseGeometry
from shapely.lib import line_merge
from shapely.ops import substring, nearest_points
from shapely.prepared import prep

from brdr.constants import BUFFER_MULTIPLICATION_FACTOR
from brdr.enums import SnapStrategy

log = logging.getLogger(__name__)


def buffer_neg_pos(geometry, buffer_value, mitre_limit=5):
    """
    Computes two buffers accordingly: one with a negative buffer value and another with
    a positive buffer value. This function can be used the check where relevant areas
    can be found (areas that will not be removed with the negative buffer), these areas
    will be positive buffered to return to their normal borders, where relevant

    Args:
        geometry (shapely.geometry.base.BaseGeometry): The input geometric object.
        buffer_value (float): The buffer distance (positive or negative).

    Returns:
        shapely.geometry.base.BaseGeometry: The union of the negative and positive
        buffers.

    Notes:
        -   The function uses the Shapely library for geometric operations.
        -   The negative buffer is applied first, followed by the positive buffer.
        -   Parameters like `quad_segs`, `join_style`, and `mitre_limit` are used for
            both buffers.

    Example:
        >>> from shapely import Point
        >>> point = Point(0, 0)
        >>> result = buffer_neg_pos(point, 1.0)
        >>> print(result)
        POLYGON EMPTY
    """
    # Implementation details
    return buffer(
        buffer(
            geometry,
            -buffer_value,
            # quad_segs=QUAD_SEGMENTS,
            join_style="mitre",
            mitre_limit=mitre_limit,
        ),
        buffer_value,
        # quad_segs=QUAD_SEGMENTS,
        join_style="mitre",
        mitre_limit=mitre_limit,
    )


def buffer_neg(geometry, buffer_value, mitre_limit=5):
    """
    Computes the negative buffer of a given geometric object.

    Args:
        geometry (shapely.geometry.base.BaseGeometry): The input geometric object.
        buffer_value (float): The negative buffer distance.

    Returns:
        shapely.geometry.base.BaseGeometry: The result of applying the negative buffer.

    Notes:
        -   The function uses the Shapely library for geometric operations.
        -   Parameters like `quad_segs`, `join_style`, and `mitre_limit` are used for
            the buffer.

    Example:
        >>> from shapely import Point
        >>> point = Point(0, 0)
        >>> result = buffer_neg(point, 1.0)
        >>> print(result)
        POLYGON EMPTY
    """
    # Implementation details
    return buffer(
        geometry,
        -buffer_value,
        # quad_segs=QUAD_SEGMENTS,
        join_style="mitre",
        mitre_limit=mitre_limit,
    )


def buffer_pos(geometry, buffer_value, mitre_limit=5):
    """
    Computes the positive buffer of a given geometric object.

    Args:
        geometry (shapely.geometry.base.BaseGeometry): The input geometric object.
        buffer_value (float): The positive buffer distance.

    Returns:
        shapely.geometry.base.BaseGeometry: The result of applying the positive buffer.

    Notes:
        -   The function uses the Shapely library for geometric operations.
        -   Parameters like `quad_segs`, `join_style`, and `mitre_limit` are used for
            the buffer.

    Example:
        >>> from shapely import Point
        >>> point = Point(0, 0)
        >>> result = buffer_pos(point, 1.0)
        >>> print(result)
        POLYGON ((-1 -1, -1 1, 1 1, 1 -1, -1 -1))
    """
    # Implementation details
    return buffer(
        geometry,
        buffer_value,
        # quad_segs=QUAD_SEGMENTS,
        join_style="mitre",
        mitre_limit=mitre_limit,
    )


def safe_union(geom_a: BaseGeometry, geom_b: BaseGeometry) -> BaseGeometry:
    """
    Safely computes the union of two geometric objects (geom_a and geom_b).

    This function handles exceptional cases where the union operation may fail due to
    non-noded intersections between LINESTRING geometries.

    Args:
        geom_a (shapely.geometry.base.BaseGeometry): The first geometric object.
        geom_b (shapely.geometry.base.BaseGeometry): The second geometric object.

    Returns:
        shapely.geometry.base.BaseGeometry: The union of geom_a and geom_b, or an empty
            Polygon if an error occurs.

    Raises:
        None

    References:
        -   Original issue: https://gis.stackexchange.com/questions/50399
        -   Buffer workaround: If the initial union fails, a small buffer is applied to
            each geometry before retrying.

    Note:
        -   The function uses the Shapely library for geometric operations.
        -   It catches exceptions related to non-noded intersections and attempts a
            buffered union.
        -   If all else fails, an empty Polygon is returned.

    Example:
        >>> from shapely import Point
        >>> point_a = Point(0, 0)
        >>> point_b = Point(1, 1)
        >>> result = safe_union(point_a, point_b)
        >>> print(result)
        POLYGON EMPTY
    """
    try:
        geom = union(geom_a, geom_b)
    except GEOSException:
        try:

            logging.warning(
                "union_error for geoms:" + geom_a.wkt + " and " + geom_b.wkt
            )
            geom = union(buffer(geom_a, 0.0000001), buffer(geom_b, 0.0000001))
        except Exception:  # noqa
            logging.error("error: empty geometry returned")
            geom = Polygon()

    return geom


def safe_equals(geom_a, geom_b):
    """
    Checks equality between two geometries with error handling.

    This function computes the equality between two Shapely geometry objects
    (`geom_a` and `geom_b`). It incorporates error handling to address potential
    exceptions that might arise due to topological inconsistencies in the
    geometries, similar to non-noded intersections between linestrings.

    Args:
        geom_a (BaseGeometry): The first Shapely geometry object.
        geom_b (BaseGeometry): The second Shapely geometry object

    Returns:
        Boolean: The equality between 2 Shapely-geometries

    Logs:
        - If a `GEOSException` occurs:
            - A warning message is logged with the WKT representations of both
                geometries.
            - The function attempts to buffer both geometries by a small value
                (0.0000001) and then perform the equality-operation.
        - If any other exception occurs:
            - An error message is logged indicating that False is returned.
    """
    # function to solve exceptional error: shapely.errors.GEOSException:
    # TopologyException: found non-noded intersection between LINESTRING
    # see: https://gis.stackexchange.com/questions/50399
    try:
        equal = equals(geom_a, geom_b)
    except GEOSException:
        logging.debug("equals_error")
        try:
            logging.warning(
                "equals_error for geoms:" + geom_a.wkt + " and " + geom_b.wkt
            )
            equal = equals(buffer(geom_a, 0.0000001), buffer(geom_b, 0.0000001))
        except Exception:  # noqa
            logging.error("equals_error: False returned")
            equal = False

    return equal


def safe_intersection(geom_a: BaseGeometry, geom_b: BaseGeometry) -> BaseGeometry:
    """
    Calculates the intersection of two geometries with error handling.

    This function attempts to compute the intersection between two Shapely geometry
    objects (`geom_a` and `geom_b`).
    It incorporates error handling to address potential exceptions that might arise
    due to topological inconsistencies
    in the geometries, such as non-noded intersections between linestrings.

    Args:
        geom_a (BaseGeometry): The first Shapely geometry object.
        geom_b (BaseGeometry): The second Shapely geometry object.

    Returns:
        BaseGeometry: The intersection geometry as a Shapely object. It might be an
        empty Polygon if an error occurs during processing.

    Logs:
        - If a `GEOSException` occurs:
            - A warning message is logged with the WKT representations of both
                geometries.
            - The function attempts to buffer both geometries by a small value
                (0.0000001) and then perform the intersection.
        - If any other exception occurs:
            - An error message is logged indicating that an empty geometry is returned.
    """
    # function to solve exceptional error: shapely.errors.GEOSException:
    # TopologyException: found non-noded intersection between LINESTRING
    # see: https://gis.stackexchange.com/questions/50399
    try:
        geom = intersection(geom_a, geom_b)
    except GEOSException:
        try:
            logging.warning(
                "intersection_error for geoms:" + geom_a.wkt + " and " + geom_b.wkt
            )
            geom = intersection(buffer(geom_a, 0.0000001), buffer(geom_b, 0.0000001))
        except Exception:  # noqa
            logging.error("error: empty geometry returned")
            geom = Polygon()

    return geom


def safe_unary_union(geometries):
    try:
        union = make_valid(unary_union(geometries))
    except:
        geometries = [make_valid(g) for g in geometries]
        union = make_valid(unary_union(geometries))

    return union


def safe_difference(geom_a, geom_b):
    """
    Calculates the difference between two geometries with error handling.

    This function computes the difference between two Shapely geometry objects
    (`geom_a` and `geom_b`). It incorporates error handling to address potential
    exceptions that might arise due to topological inconsistencies in the
    geometries, similar to non-noded intersections between linestrings.

    Args:
        geom_a (BaseGeometry): The first Shapely geometry object.
        geom_b (BaseGeometry): The second Shapely geometry object to be subtracted from
        the first.

    Returns:
        BaseGeometry: The difference geometry as a Shapely object. It might be an empty
        Polygon if an error occurs during processing.

    Logs:
        - If a `GEOSException` occurs:
            - A warning message is logged with the WKT representations of both
                geometries.
            - The function attempts to buffer both geometries by a small value
                (0.0000001) and then perform the difference operation.
        - If any other exception occurs:
            - An error message is logged indicating that an empty geometry is returned.
    """
    # function to solve exceptional error: shapely.errors.GEOSException:
    # TopologyException: found non-noded intersection between LINESTRING
    # see: https://gis.stackexchange.com/questions/50399
    try:
        geom = difference(geom_a, geom_b)
    except GEOSException:
        logging.debug("difference_error")
        try:
            logging.warning(
                "difference_error for geoms:" + geom_a.wkt + " and " + geom_b.wkt
            )
            geom = difference(buffer(geom_a, 0.0000001), buffer(geom_b, 0.0000001))
        except Exception:  # noqa
            logging.error("error: empty geometry returned")
            geom = Polygon()

    return geom


def safe_symmetric_difference(geom_a, geom_b):
    """
    Calculates the symmetrical difference between two geometries with error handling.

    This function computes the symmetrical difference between two Shapely geometry
    objects (`geom_a` and `geom_b`). It incorporates error handling to address
    potential exceptions that might arise due to topological inconsistencies in
    the geometries, similar to non-noded intersections between linestrings.

    Args:
        geom_a (BaseGeometry): The first Shapely geometry object.
        geom_b (BaseGeometry): The second Shapely geometry object to be subtracted
        from the first.

    Returns:
        BaseGeometry: The symmetrical difference geometry as a Shapely object. It
        might be an empty Polygon if an error occurs during processing.

    Logs:
        - If a `GEOSException` occurs:
            - A warning message is logged with the WKT representations of both
              geometries.
            - The function attempts to buffer both geometries by a small value (
              0.0000001) and then perform the difference operation.
        - If any other exception occurs:
            - An error message is logged indicating that an empty geometry is returned.
    """
    # function to solve exceptional error: shapely.errors.GEOSException:
    # TopologyException: found non-noded intersection between LINESTRING
    # see: https://gis.stackexchange.com/questions/50399
    try:
        geom = symmetric_difference(geom_a, geom_b)
    except GEOSException:
        try:
            logging.warning(
                "symmetric_difference_error for geoms:"
                + geom_a.wkt
                + " and "
                + geom_b.wkt
            )
            geom = symmetric_difference(
                buffer(geom_a, 0.0000001), buffer(geom_b, 0.0000001)
            )
        except Exception:  # noqa
            logging.error("error: empty geometry returned")
            geom = Polygon()

    return geom


def _grid_bounds(geom: BaseGeometry, delta: float):
    """
    Divides a geometric area (specified by `geom`) into a grid of rectangular
    partitions.

    Args:
        geom (BaseGeometry): The geometric object representing the area to be
            partitioned.
        delta (float): The desired distance between partitions.

    Returns:
        list: A list of Polygon objects representing the divided partitions.
    """
    if geom is None or is_empty(geom):
        return geom
    bounds = geom.bounds
    min_x, min_y, max_x, max_y = bounds
    nx = int((max_x - min_x) / delta)
    ny = int((max_y - min_y) / delta)
    if nx < 2 and ny < 2:
        return [geom.envelope]
    if nx < 2:
        nx = 2
    if ny < 2:
        ny = 2
    gx, gy = np.linspace(min_x, max_x, nx + 1), np.linspace(min_y, max_y, ny + 1)
    grid = []
    for i in range(len(gx) - 1):
        for j in range(len(gy) - 1):
            poly_ij = Polygon(
                [
                    [gx[i], gy[j]],
                    [gx[i], gy[j + 1]],
                    [gx[i + 1], gy[j + 1]],
                    [gx[i + 1], gy[j]],
                ]
            )
            grid.append(poly_ij)
    return grid


def geom_from_wkt(wkt_string):
    """
    Converts a WellKnownText (WKT) into a shapely-geometry
    Args:
        wkt_string:  WellKnownText (WKT)

    Returns: Shapely geometry

    """
    return from_wkt(wkt_string)


def geom_to_wkt(shapely_geometry):
    """
    Converts a shapely-geometry into WellKnownText (WKT)
    Args:
        shapely_geometry:

    Returns: WellKnownText (WKT) - string

    """
    return to_wkt(shapely_geometry)


def snap_geometry_to_reference(
    geometry,
    reference,
    snap_strategy=SnapStrategy.NO_PREFERENCE,
    max_segment_length=-1,
    tolerance=1,
):
    if geometry is None or geometry.is_empty or reference is None or reference.is_empty:
        return geometry
    if max_segment_length > 0:
        geometry = segmentize(geometry, max_segment_length=max_segment_length)
    geometry = to_multi(geometry, geomtype=None)
    if geometry.geom_type == "MultiPoint":
        result = _snap_point_to_reference(
            geometry,
            reference=reference,
            snap_strategy=snap_strategy,
            max_segment_length=max_segment_length,
            tolerance=tolerance,
        )
    elif geometry.geom_type == "MultiLineString":
        result = _snap_line_to_reference(
            geometry,
            reference=reference,
            snap_strategy=snap_strategy,
            max_segment_length=max_segment_length,
            tolerance=tolerance,
        )
    elif geometry.geom_type == "MultiPolygon":
        result = _snap_polygon_to_reference(
            geometry,
            reference=reference,
            snap_strategy=snap_strategy,
            max_segment_length=max_segment_length,
            tolerance=tolerance,
        )

    elif geometry.geom_type == "GeometryCollection":
        results = []
        for geom in geometry.geoms:
            result = snap_geometry_to_reference(
                geom,
                reference=reference,
                snap_strategy=snap_strategy,
                max_segment_length=max_segment_length,
                tolerance=tolerance,
            )
            results.append(result)
        # result = GeometryCollection(results)
        result = safe_unary_union(results)

    else:
        raise NotImplementedError(
            f"snapping for this type of geometry is not implemented: {str(geometry.geom_type)}"
        )
    return result


def _snap_point_to_reference(
    geometry,
    reference,
    snap_strategy=SnapStrategy.NO_PREFERENCE,
    max_segment_length=-1,
    tolerance=1,
):
    if geometry is None or geometry.is_empty or reference is None or reference.is_empty:
        return geometry
    if max_segment_length > 0:
        geometry = segmentize(geometry, max_segment_length=max_segment_length)

    if geometry.geom_type == "Point":
        geometry = MultiPoint([geometry])

    ref_border, ref_borders, ref_coords = _get_ref_objects(reference)

    if len(ref_coords) == 0:
        snap_strategy = SnapStrategy.NO_PREFERENCE

    points = []
    for geom in geometry.geoms:
        coords = list(geom.coords)
        coordinates = []
        for idx, coord in enumerate(coords):
            p = Point(coords[idx])
            p_snapped, bool_snapped, ref_vertices = _get_snapped_point(
                p, ref_border, ref_coords, snap_strategy, tolerance
            )
            coordinates.append(p_snapped.coords[0])

        # convert coordinates back to a point
        point = make_valid(Point(coordinates))
        points.append((point))
    result = safe_unary_union(points)
    return result


def _snap_line_to_reference(
    geometry,
    reference,
    snap_strategy=SnapStrategy.NO_PREFERENCE,
    max_segment_length=-1,
    tolerance=1,
):
    # return snap(geometry,reference,tolerance)
    if geometry is None or geometry.is_empty or reference is None or reference.is_empty:
        return geometry
    if max_segment_length > 0:
        geometry = segmentize(geometry, max_segment_length=max_segment_length)

    if geometry.geom_type == "LineString":
        geometry = MultiLineString([geometry])

    ref_border, ref_borders, ref_coords = _get_ref_objects(reference)
    if len(ref_coords) == 0:
        snap_strategy = SnapStrategy.NO_PREFERENCE

    lines = []
    for geom in geometry.geoms:
        coords = list(geom.coords)
        coordinates = _get_snapped_coordinates(
            coords, ref_border, ref_borders, ref_coords, snap_strategy, tolerance
        )

        # convert coordinates back to a line
        linestring = make_valid(LineString(coordinates))
        lines.append((linestring))
    result = safe_unary_union(lines)
    result = longest_linestring_from_multilinestring(
        result
    )  # Solves overshoots in multilinestrings to return a linestring
    return result


def _snap_polygon_to_reference(
    geometry,
    reference,
    snap_strategy=SnapStrategy.PREFER_VERTICES,
    max_segment_length=-1,
    tolerance=1,
    correction_distance=0.01,
):
    if geometry is None or geometry.is_empty or reference is None or reference.is_empty:
        return geometry
    if max_segment_length > 0:
        geometry = segmentize(geometry, max_segment_length=max_segment_length)
    geometry = to_multi(geometry, geomtype="Polygon")

    ref_border, ref_borders, ref_coords = _get_ref_objects(reference)
    if len(ref_coords) == 0:
        snap_strategy = SnapStrategy.NO_PREFERENCE
    polygons = []
    for geom in geometry.geoms:
        coords = list(geom.exterior.coords)
        coordinates = _get_snapped_coordinates(
            coords, ref_border, ref_borders, ref_coords, snap_strategy, tolerance
        )
        # convert coordinates back to a polygon
        polygon = make_valid(Polygon(coordinates))
        polygons.append((polygon))
    return buffer_neg_pos(safe_unary_union(polygons), correction_distance)


def _get_ref_objects(reference):
    # reference = safe_intersection(
    #     reference, buffer_pos(geometry, tolerance * BUFFER_MULTIPLICATION_FACTOR)
    # )
    ref_coords = list(get_coords_from_geometry(reference))
    reference_list = []
    if reference.geom_type == "GeometryCollection":
        for r in reference.geoms:
            reference_list.append(to_multi(r, geomtype=None))
    else:
        reference_list.append(to_multi(reference, geomtype=None))
    ref_borders = []
    for r in reference_list:
        for g in r.geoms:
            if g.geom_type == "Polygon":
                ref_borders.append(g.exterior)
                ref_borders.extend([interior for interior in g.interiors])
            elif g.geom_type == "LineString":
                ref_borders.append(g)
            elif g.geom_type == "Point":
                ref_borders.append(g)  # point
    ref_border = safe_unary_union(ref_borders)
    return ref_border, ref_borders, ref_coords


def _get_snapped_coordinates(
    coords, ref_border, ref_borders, ref_coords, snap_strategy, tolerance
):
    coordinates = []
    for idx, coord in enumerate(coords):  # for each vertex in the first line
        if idx == 0:
            continue
        p_start = Point(coords[idx - 1])
        p_end = Point(coords[idx])

        p_start_snapped, bool_start_snapped, ref_vertices_start = _get_snapped_point(
            p_start, ref_border, ref_coords, snap_strategy, tolerance
        )
        p_end_snapped, bool_end_snapped, ref_vertices_end = _get_snapped_point(
            p_end, ref_border, ref_coords, snap_strategy, tolerance
        )

        coordinates.append(p_start_snapped.coords[0])

        if not bool_start_snapped and not bool_end_snapped:
            coordinates.append(p_end_snapped.coords[0])
            continue
        elif bool_start_snapped and bool_end_snapped:
            coordinates.extend(
                _get_sublinestring_coordinates(
                    p_end,
                    p_end_snapped,
                    p_start,
                    p_start_snapped,
                    ref_borders,
                    tolerance,
                )
            )
            coordinates.append(p_end_snapped.coords[0])
        elif bool_start_snapped != bool_end_snapped:
            # determine_p_mid
            line = LineString([p_start, p_end])
            # print ("idx:" + str(idx))
            # print (line.wkt)
            if ref_vertices_start or ref_vertices_end:
                ref_buffered = buffer_pos(MultiPoint(ref_coords), tolerance)
            else:
                ref_buffered = buffer_pos(ref_border, tolerance)
            intersected_line = safe_intersection(line, ref_buffered)
            if intersected_line.is_empty or intersected_line.geom_type != "LineString":
                coordinates.append(p_end_snapped.coords[0])
                continue
            intersected_line_boundary_points = intersected_line.boundary.geoms
            first = intersected_line_boundary_points[0]
            last = intersected_line_boundary_points[-1]
            if first == p_start:
                p_mid = last
                p_mid_snapped, bool_mid_snapped, ref_vertices_mid = _get_snapped_point(
                    p_mid, ref_border, ref_coords, snap_strategy, tolerance
                )

                coordinates.extend(
                    _get_sublinestring_coordinates(
                        p_mid,
                        p_mid_snapped,
                        p_start,
                        p_start_snapped,
                        ref_borders,
                        tolerance,
                    )
                )
                coordinates.append(p_mid_snapped.coords[0])
                coordinates.append(p_end_snapped.coords[0])
            elif last == p_end:
                p_mid = first
                p_mid_snapped, bool_mid_snapped, ref_vertices_mid = _get_snapped_point(
                    p_mid, ref_border, ref_coords, snap_strategy, tolerance
                )
                coordinates.append(p_mid_snapped.coords[0])
                coordinates.extend(
                    _get_sublinestring_coordinates(
                        p_end,
                        p_end_snapped,
                        p_mid,
                        p_mid_snapped,
                        ref_borders,
                        tolerance,
                    )
                )
                coordinates.append(p_end_snapped.coords[0])
    return coordinates


def _get_sublinestring_coordinates(
    p_end, p_end_snapped, p_start, p_start_snapped, ref_borders, tolerance
):
    coordinates = []
    if p_start_snapped == p_end_snapped:
        return coordinates
    reference_border, distance = closest_line(ref_borders, p_end)
    if reference_border is None or reference_border.is_empty:
        return coordinates
    distance_start_end = p_start.distance(p_end)
    line_substring = _get_line_substring(
        reference_border, p_start_snapped, p_end_snapped, distance_start_end
    )
    if line_substring is None or line_substring.is_empty:
        return coordinates
    # TODO what if line_substring is a multilinestring? error detected in brdrQ
    for p in line_substring.coords:
        point = Point(p)
        if (point.distance(p_start) + point.distance(p_end)) / 2 <= tolerance:
            coordinates.append(p)
    return coordinates


def geometric_equality(geom_a, geom_b, correction_distance, mitre_limit):
    return buffer_neg(
        safe_symmetric_difference(geom_a, geom_b),
        correction_distance,
        mitre_limit=mitre_limit,
    ).is_empty


def _get_line_substring(
    reference_border, p_start_snapped, p_end_snapped, distance_start_end
):
    if (
        reference_border is None
        or reference_border.is_empty
        or reference_border.geom_type
        not in ["LinearRing", "LineString", "MultiLineString"]
    ):
        return reference_border
    start_fraction = reference_border.project(p_start_snapped, normalized=True)
    end_fraction = reference_border.project(p_end_snapped, normalized=True)
    line_substring = substring(
        reference_border,
        start_dist=start_fraction,
        end_dist=end_fraction,
        normalized=True,
    )
    if line_substring.length > 1.5 * distance_start_end:

        if start_fraction < end_fraction:
            subline_a = substring(
                reference_border, start_dist=start_fraction, end_dist=0, normalized=True
            )
            subline_b = substring(
                reference_border, start_dist=100, end_dist=end_fraction, normalized=True
            )
        else:
            subline_a = substring(
                reference_border,
                start_dist=start_fraction,
                end_dist=100,
                normalized=True,
            )
            subline_b = substring(
                reference_border, start_dist=0, end_dist=end_fraction, normalized=True
            )

        sublines = [s for s in [subline_a, subline_b] if s.geom_type == "LineString"]
        line_substring_2 = line_merge(MultiLineString(sublines))
        if (
            line_substring_2.length < line_substring.length
            and line_substring_2.geom_type == ["LineString"]
        ):
            line_substring = line_substring_2
        if line_substring.length > 3 * distance_start_end:
            line_substring = LineString([p_start_snapped, p_end_snapped])
    return line_substring


def _snapped_point_by_snapstrategy(
    p, p_nearest, p_nearest_vertices, snap_strategy, tolerance
):
    p_snapped = p
    snapped = False
    ref_vertices = False
    if p_nearest is None:
        return p_snapped, snapped, ref_vertices
    if snap_strategy == SnapStrategy.NO_PREFERENCE:
        if p.distance(p_nearest) <= tolerance:
            p_snapped = p_nearest
            snapped = True
        else:
            p_snapped = p
    elif snap_strategy == SnapStrategy.ONLY_VERTICES:
        if (
            p_nearest_vertices is not None
            and p.distance(p_nearest_vertices) <= tolerance
        ):
            p_snapped = p_nearest_vertices
            snapped = True
            ref_vertices = True
        else:
            p_snapped = p
    elif snap_strategy == SnapStrategy.PREFER_VERTICES:
        if (
            p_nearest_vertices is not None
            and p.distance(p_nearest_vertices) <= tolerance
        ):
            p_snapped = p_nearest_vertices
            snapped = True
            ref_vertices = True
        elif p.distance(p_nearest) <= tolerance:
            p_snapped = p_nearest
            snapped = True
        else:
            p_snapped = p
    return p_snapped, snapped, ref_vertices


def closest_line(lines, point):
    # get distances
    distance_list = [line.distance(point) for line in lines]
    if len(distance_list) == 0:
        return None, -1
    shortest_distance = min(distance_list)  # find the line closest to the point
    return (
        lines[distance_list.index(shortest_distance)],  # return the closest line
        shortest_distance,
    )  # return the distance to that line


def create_donut(geometry, distance):
    if distance == 0:
        return geometry
    inner_geometry = buffer_neg(geometry, distance)
    return safe_difference(geometry, inner_geometry)


def features_by_geometric_operation(
    list_input_geometries, list_input_ids, list_geometries, predicate="intersects"
):
    tree = STRtree(list_input_geometries)
    id_items = np.array(list_input_ids)
    arr_indices = tree.query(list_geometries, predicate=predicate)
    resulting_id_items = list(set(id_items.take(arr_indices[1])))
    # resulting_id_items = [str(element) for element in resulting_id_items]
    return resulting_id_items


def get_partitions(geom, delta):
    """
    Filters a computed grid of partitions (generated by `_grid_bounds`) based on
    intersection with a geometric object (`geom`).

    Args:
        geom (BaseGeometry): The geometric object to check for intersection
            with partitions.
        delta (float): The distance between partitions (same value used in
            `_grid_bounds`).

    Returns:
        list: A filtered list of Polygon objects representing the partitions
            overlapping the original geometric object.
    """
    # TODO: partitioning results in multiple squares, this can be improved by
    #  partitioning with a quadtree with rectangles?
    #  https://www.fundza.com/algorithmic/quadtree/index.html
    prepared_geom = prep(geom)
    partitions = _grid_bounds(geom, delta)
    filtered_grid = list(filter(prepared_geom.intersects, partitions))
    return filtered_grid


def get_shape_index(area, perimeter):
    if area > 0 and perimeter > 0:
        return 4 * pi * (area / (perimeter**2))
    else:
        return -1


def fill_and_remove_gaps(input_geometry, buffer_value):
    cleaned_geometry = input_geometry
    ix_part = 1
    for part in get_parts(input_geometry):
        exterior_ring = get_exterior_ring(part)
        exterior_polygon = polygons([exterior_ring])[0]
        empty_buffered_exterior_polygon = buffer_neg(
            exterior_polygon, buffer_value
        ).is_empty
        if (
            ix_part > 1
            and empty_buffered_exterior_polygon
            and not exterior_polygon.is_empty
        ):
            cleaned_geometry = safe_difference(cleaned_geometry, exterior_polygon)
        num_interior_rings = get_num_interior_rings(part)
        if num_interior_rings > 0:
            ix = 0
            while ix < num_interior_rings:
                interior_ring = get_interior_ring(part, ix)
                interior_polygon = polygons([interior_ring])[0]

                empty_buffered_interior_ring = buffer_neg(
                    interior_polygon, buffer_value
                ).is_empty
                if empty_buffered_interior_ring:
                    cleaned_geometry = safe_union(cleaned_geometry, interior_polygon)
                ix = ix + 1
        ix_part = ix_part + 1
    return cleaned_geometry


def get_bbox(geometry):
    """
    Get the BBOX (string) of a shapely geometry
    """
    return str(geometry.bounds).strip("()")


def geojson_to_multi(geojson):
    """
    Transforms a geojson: Checks if there are single-geometry-features and transforms them into Multi-geometries, so all objects are of type 'Multi' (or null-geometry).
    It is important that geometry-type is consitent (f.e. in QGIS) to show and style the geojson-layer
    """

    if geojson is None or "features" not in geojson or geojson["features"] is None:
        return geojson
    for f in geojson["features"]:
        if f["geometry"] is None:
            continue
        if f["geometry"]["type"] == "Polygon":
            f["geometry"] = {
                "type": "MultiPolygon",
                "coordinates": [f["geometry"]["coordinates"]],
            }
        elif f["geometry"]["type"] == "LineString":
            f["geometry"] = {
                "type": "MultiLineString",
                "coordinates": [f["geometry"]["coordinates"]],
            }
        elif f["geometry"]["type"] == "Point":
            f["geometry"] = {
                "type": "MultiPoint",
                "coordinates": [f["geometry"]["coordinates"]],
            }
    return geojson


def to_multi(geometry, geomtype=None):
    """
    Converts a Shapely geometry to its corresponding multi-variant.

    Parameters:
    geometry (shapely.geometry.base.BaseGeometry): The input geometry to be converted.
    geomtype (str, optional): The type of geometry to extract and convert.
                              Possible values are 'Point', 'LineString', 'Polygon', or None.
                              If None, all geometries are converted to their multi-variant.

    Returns:
    shapely.geometry.base.BaseGeometry: The converted multi-geometry or an empty multi-geometry of the specified type.

    Raises:
    TypeError: If an unknown geometry type is provided.

    Examples:
    >>> from shapely.geometry import Point, LineString, Polygon, GeometryCollection
    >>> point = Point(1, 1)
    >>> line = LineString([(0, 0), (1, 1)])
    >>> polygon = Polygon([(0, 0), (1, 1), (1, 0)])
    >>> collection = GeometryCollection([point, line, polygon])
    >>> to_multi(point, 'Point')
    <shapely.geometry.multipoint.MultiPoint object at 0x...>
    >>> to_multi(line, 'LineString')
    <shapely.geometry.multilinestring.MultiLineString object at 0x...>
    >>> to_multi(polygon, 'Polygon')
    <shapely.geometry.multipolygon.MultiPolygon object at 0x...>
    >>> to_multi(collection, 'Point')
    <shapely.geometry.multipoint.MultiPoint object at 0x...>
    >>> to_multi(collection)
    <shapely.geometry.collection.GeometryCollection object at 0x...>
    """
    if geometry is None:
        return geometry
    if geomtype in ["Point", "MultiPoint"]:
        if isinstance(geometry, Point):
            return MultiPoint([geometry])
        elif isinstance(geometry, MultiPoint):
            return geometry
        elif isinstance(geometry, GeometryCollection):
            points = [geom for geom in geometry.geoms if isinstance(geom, Point)]
            return MultiPoint(points)
        else:
            return MultiPoint()
    elif geomtype in ["LineString", "MultiLineString"]:
        if isinstance(geometry, LineString):
            return MultiLineString([geometry])
        elif isinstance(geometry, MultiLineString):
            return geometry
        elif isinstance(geometry, GeometryCollection):
            lines = [geom for geom in geometry.geoms if isinstance(geom, LineString)]
            return MultiLineString(lines)
        else:
            return MultiLineString()
    elif geomtype in ["Polygon", "MultiPolygon"]:
        if isinstance(geometry, Polygon):
            return MultiPolygon([geometry])
        elif isinstance(geometry, MultiPolygon):
            return geometry
        elif isinstance(geometry, GeometryCollection):
            polygons = [geom for geom in geometry.geoms if isinstance(geom, Polygon)]
            return MultiPolygon(polygons)
        else:
            return MultiPolygon()
    elif geomtype is None:
        if isinstance(geometry, Point):
            if not geometry.is_empty:
                return MultiPoint([geometry])
            else:
                return MultiPoint()
        elif isinstance(geometry, LineString):
            if not geometry.is_empty:
                return MultiLineString([geometry])
            else:
                return MultiLineString()
        elif isinstance(geometry, Polygon):
            if not geometry.is_empty:
                return MultiPolygon([geometry])
            else:
                return MultiPolygon()
        elif (
            isinstance(geometry, MultiPoint)
            or isinstance(geometry, MultiLineString)
            or isinstance(geometry, MultiPolygon)
        ):
            return geometry  # Het is al een multi-variant
        elif isinstance(geometry, GeometryCollection):
            multi_geoms = []
            for geom in geometry.geoms:
                if isinstance(geom, Point):
                    multi_geoms.append(MultiPoint([geom]))
                elif isinstance(geom, LineString):
                    multi_geoms.append(MultiLineString([geom]))
                elif isinstance(geom, Polygon):
                    multi_geoms.append(MultiPolygon([geom]))
                elif (
                    isinstance(geom, MultiPoint)
                    or isinstance(geom, MultiLineString)
                    or isinstance(geom, MultiPolygon)
                ):
                    multi_geoms.append(geom)
            return GeometryCollection(multi_geoms)
        else:
            raise TypeError("Geometry type not supported: {}".format(type(geometry)))
    else:
        raise TypeError("Geometry type not supported: {}".format(type(geometry)))


def get_coords_from_geometry(geometry):
    coords = set()
    if geometry is None or geometry.is_empty:
        return coords
    elif isinstance(geometry, Point):
        coords.update(geometry.coords)
    elif isinstance(geometry, MultiPoint):
        for pt in geometry.geoms:
            coords.update(get_coords_from_geometry(pt))
    elif isinstance(geometry, LineString):
        coords.update(geometry.coords)
    elif isinstance(geometry, MultiLineString):
        for line in geometry.geoms:
            coords.update(get_coords_from_geometry(line))
    elif isinstance(geometry, Polygon):
        coords.update(geometry.exterior.coords)
        for linearring in geometry.interiors:
            coords.update(linearring.coords)
    elif isinstance(geometry, MultiPolygon):
        for polygon in geometry.geoms:
            coords.update(get_coords_from_geometry(polygon))
    elif isinstance(geometry, GeometryCollection):
        for geom in geometry.geoms:
            coords.update(get_coords_from_geometry(geom))
    return coords


def get_geoms_from_geometry(geometry):
    geoms = set()
    if geometry is None or geometry.is_empty:
        return geoms
    elif isinstance(geometry, (Point, LineString, Polygon)):
        geoms.update([geometry])
    elif isinstance(geometry, (MultiPoint, MultiLineString, MultiPolygon)):
        geoms.update(geometry.geoms)
    elif isinstance(geometry, GeometryCollection):
        for geom in geometry.geoms:
            geoms.update(get_geoms_from_geometry(geom))
    return geoms


def snap_multilinestring_endpoints(multilinestring, tolerance):
    """
    Snapt de begin- en eindpunten van elke LineString in een MultiLineString
    aan nabijgelegen lijnen binnen een gegeven tolerantie.

    Parameters:
    - multilinestring: shapely.geometry.MultiLineString
    - tolerance: float, afstandstolerantie voor snappen

    Returns:
    - Een nieuwe MultiLineString met gesnapte eindpunten
    """
    if isinstance(multilinestring, LineString):
        return MultiLineString([multilinestring])
    lines = list(multilinestring.geoms)
    snapped_lines = []

    for i, line in enumerate(lines):
        start = Point(line.coords[0])
        end = Point(line.coords[-1])

        # Verzamel alle andere lijnen
        other_lines = [l for j, l in enumerate(lines) if j != i]

        # Snap begin- en eindpunt indien binnen tolerantie
        for other in other_lines:
            if start.distance(other) <= tolerance:
                start = snap(start, other, tolerance)
            if end.distance(other) <= tolerance:
                end = snap(end, other, tolerance)

        # Herbouw de lijn met eventueel gesnapte punten
        new_coords = [start.coords[0]] + list(line.coords[1:-1]) + [end.coords[0]]
        snapped_lines.append(LineString(new_coords))

    return MultiLineString(snapped_lines)


def get_connection_lines_to_nearest(multilinestring):
    # Extract all endpoints
    endpoints = []
    for line in multilinestring.geoms:
        coords = list(line.coords)
        endpoints.append(Point(coords[0]))
        endpoints.append(Point(coords[-1]))

    # Count occurrences of each endpoint
    endpoint_counts = Counter((pt.x, pt.y) for pt in endpoints)

    # Identify loose endpoints (those that appear only once)
    loose_endpoints = [Point(xy) for xy, count in endpoint_counts.items() if count == 1]

    # Keep track of which points have been connected
    used = set()
    connection_lines = []

    while loose_endpoints:
        pt = loose_endpoints.pop(0)
        if (pt.x, pt.y) in used:
            continue
        # Find the closest other loose endpoint
        closest_pt = min(
            (other for other in loose_endpoints if (other.x, other.y) not in used),
            key=lambda p: pt.distance(p),
            default=None,
        )
        if closest_pt:
            connection_lines.append(LineString([pt, closest_pt]))
            used.add((pt.x, pt.y))
            used.add((closest_pt.x, closest_pt.y))
            loose_endpoints.remove(closest_pt)

    # Combine original lines with connection lines
    # all_lines = list(multilinestring.geoms) + connection_lines
    # return MultiLineString(all_lines)
    return connection_lines


def shortest_connections_between_geometries(geometry):
    """
    Voor elk element in een GeometryCollection, bepaal de kortste verbindingslijn
    naar het dichtstbijzijnde andere element.

    Parameters:
    - geometry_collection: shapely.geometry.GeometryCollection

    Returns:
    - List van shapely.geometry.LineString objecten die de kortste verbindingen representeren
    """
    connection_lines = []
    geometry = safe_unary_union(geometry)
    geometry = to_multi(geometry)

    if geometry is None or geometry.is_empty or len(geometry.geoms) <= 1:
        return connection_lines

    geometries = list(geometry.geoms)

    for i, geom in enumerate(geometries):
        other_geometries = list(geometry.geoms)
        del other_geometries[i]
        connection_line = shortest_line(geom, safe_unary_union(other_geometries))
        connection_lines.append(connection_line)

    return make_linestrings_unique(connection_lines)


def make_linestrings_unique(lines):
    """
    Verwijder dubbele LineStrings uit een lijst, waarbij lijnen met omgekeerde richting
    als gelijk worden beschouwd.

    Parameters:
    - lines: lijst van shapely.geometry.LineString objecten

    Returns:
    - lijst van unieke LineStrings
    """
    seen = set()
    unique_lines = []

    for line in lines:
        coords = tuple(line.coords)
        # Sorteer de coördinaten zodat richting niet uitmaakt
        key = tuple(sorted([coords[0], coords[-1]]))
        if key not in seen:
            seen.add(key)
            unique_lines.append(line)

    return unique_lines


def fill_gaps_in_multilinestring(multilinestring, tolerance):
    """
    Vul kleine gaps in een MultiLineString op door verbindingslijnen toe te voegen
    tussen eindpunten die binnen een bepaalde tolerantie van elkaar liggen.

    Parameters:
    - multilinestring: shapely.geometry.MultiLineString
    - tolerance: float, maximale afstand tussen eindpunten om een verbindingslijn toe te voegen

    Returns:
    - Een nieuwe MultiLineString met toegevoegde verbindingslijnen
    """

    if isinstance(multilinestring, LineString):
        return MultiLineString([multilinestring])
    lines = list(multilinestring.geoms)
    endpoints = []

    # Verzamel alle begin- en eindpunten
    for line in lines:
        endpoints.append(Point(line.coords[0]))
        endpoints.append(Point(line.coords[-1]))

    # Zoek paren van eindpunten die dicht bij elkaar liggen
    new_lines = []
    used_pairs = set()
    for p1, p2 in combinations(endpoints, 2):
        if p1.equals(p2) or (p1, p2) in used_pairs or (p2, p1) in used_pairs:
            continue
        if p1.distance(p2) <= tolerance:
            new_lines.append(LineString([p1, p2]))
            used_pairs.add((p1, p2))

    # Voeg originele lijnen en nieuwe verbindingslijnen samen
    all_lines = lines + new_lines
    return line_merge(all_lines)


def nearest_node(point, nodes):
    """
    vind de dichtstbijzijnde knopen bij start- en eindpunt
    :param point:
    :param nodes:
    :return:
    """
    return min(nodes, key=lambda n: Point(n).distance(point))


def find_best_circle_path(directed_graph, geom_to_process):
    cycles = list(nx.simple_cycles(directed_graph))
    min_dist = inf
    best_cycle_line = None
    for cycle in cycles:
        cycle_coords = cycle + [cycle[0]]
        cycle_line = LineString(cycle_coords)
        dist = total_vertex_distance(cycle_line, geom_to_process)
        if dist < min_dist:
            min_dist = dist
            best_cycle_line = cycle_line
    return best_cycle_line


def find_circle_path(directed_graph):
    # Vind alle eenvoudige cycli
    cycles = list(nx.simple_cycles(directed_graph))

    # Bepaal de langste cyclus op basis van gewichten
    def cycle_weight(cycle):
        weight = 0
        for i in range(len(cycle)):
            u = cycle[i]
            v = cycle[(i + 1) % len(cycle)]
            if directed_graph.has_edge(u, v):
                weight += directed_graph[u][v]["weight"]
            else:
                return -1  # ongeldig pad
        return weight

    # Selecteer de langste cyclus
    longest_cycle = max(cycles, key=cycle_weight)
    # Zet de cyclus om naar een gesloten LineString
    longest_cycle_coords = longest_cycle + [longest_cycle[0]]
    longest_cycle_linestring = LineString(longest_cycle_coords)
    if not isinstance(safe_unary_union(Polygon(longest_cycle_linestring)), Polygon):
        longest_cycle_linestring = None
    return longest_cycle_linestring


def longest_linestring_from_multilinestring(multilinestring):
    """searches for a linestring in a multilinestring (with no loops)"""

    if multilinestring is None or multilinestring.is_empty:
        return GeometryCollection()

    if isinstance(multilinestring, LineString):
        return multilinestring

    if not isinstance(multilinestring, MultiLineString):
        logging.warning(
            "Multilinestring expected. Other type detected; empty geometry returned"
        )
        return GeometryCollection()

    # Create a graph from the MultiLineString
    G = graph_from_multilinestring(multilinestring)

    # Find all simple paths and keep the longest one
    longest_path = []
    max_length = 0

    for source in G.nodes:
        for target in G.nodes:
            if source != target:
                try:
                    path = nx.shortest_path(
                        G, source=source, target=target, weight="weight"
                    )
                    length = sum(
                        LineString([path[i], path[i + 1]]).length
                        for i in range(len(path) - 1)
                    )
                    if length > max_length:
                        max_length = length
                        longest_path = path
                except nx.NetworkXNoPath:
                    continue

    # Convert the longest path to a LineString
    return LineString(longest_path)


def find_longest_path_in_network(
    geom_to_process, multilinestring, snap_strategy, relevant_distance
):
    """
    Bepaal het langste pad tussen twee punten in een MultiLineString-geometrie.

    Parameters:
    - multilinestring: shapely.geometry.MultiLineString
    - start_point: shapely.geometry.Point
    - end_point: shapely.geometry.Point

    Returns:
    - shapely.geometry.LineString van het langste pad tussen de twee punten
    """
    if multilinestring is None or multilinestring.is_empty:
        return GeometryCollection()

    if isinstance(multilinestring, LineString):
        return multilinestring

    start_point = Point(geom_to_process.coords[0])
    end_point = Point(geom_to_process.coords[-1])

    G = graph_from_multilinestring(multilinestring)

    start_node = nearest_node(start_point, G.nodes)
    end_node = nearest_node(end_point, G.nodes)
    if start_node == end_node:
        return find_circle_path(G.to_directed())
    if not snap_strategy is None and snap_strategy != SnapStrategy.NO_PREFERENCE:
        # when SnapStrategy = PREFER_VERTICES or ONLY_VERTICES
        start_node = get_vertex_node(G, relevant_distance, start_node, start_point)
        end_node = get_vertex_node(G, relevant_distance, end_node, end_point)

    all_paths = nx.all_simple_paths(G, source=start_node, target=end_node)

    max_length = 0
    longest_path = []
    # TODO; when having long lines, this loop becomes very slow due to amount of possibities to check
    for path in all_paths:
        length = sum(G[path[i]][path[i + 1]]["weight"] for i in range(len(path) - 1))
        if length > max_length:
            max_length = length
            longest_path = path

    if longest_path:
        return LineString(longest_path)
    else:
        return None


def graph_from_multilinestring(multilinestring):
    if not isinstance(multilinestring, MultiLineString):
        raise TypeError("multilinstring expected")
    G = nx.Graph()
    for line in multilinestring.geoms:
        coords = list(line.coords)
        for i in range(len(coords) - 1):
            p1 = tuple(coords[i])
            p2 = tuple(coords[i + 1])
            segment = LineString([p1, p2])
            G.add_edge(p1, p2, weight=segment.length, geometry=segment)
        G, added_edges = connect_components_greedy(G)
    return G


def euclidean_distance(p1, p2):
    return Point(p1).distance(Point(p2))


def connect_components_greedy(G):
    added_edges = []

    while not nx.is_connected(G):
        components = list(nx.connected_components(G))
        min_dist = float("inf")
        best_pair = None

        # Compare nodes from different components
        for i in range(len(components)):
            for j in range(i + 1, len(components)):
                comp1 = components[i]
                comp2 = components[j]
                for u in comp1:
                    for v in comp2:
                        dist = euclidean_distance(u, v)
                        if dist < min_dist:
                            min_dist = dist
                            best_pair = (u, v)

        # Add the shortest edge found
        if best_pair:
            G.add_edge(best_pair[0], best_pair[1], weight=min_dist)
            added_edges.append((best_pair[0], best_pair[1], min_dist))

    return G, added_edges


def connect_disconnected_networks(G):
    """
    Verbindt de losse componenten in een NetworkX-graaf met de kortste mogelijke verbindingen.
    Retourneert een nieuwe graaf waarin alle componenten verbonden zijn.
    """
    # Maak een kopie van de graaf
    G_connected = G.copy()

    # Vind alle verbonden componenten
    components = list(nx.connected_components(G_connected))

    # Stop als de graaf al verbonden is
    if len(components) <= 1:
        return G_connected

    # Bepaal representatieve knopen (alle eindpunten) per component
    component_endpoints = []
    for comp in components:
        endpoints = [node for node in comp if G_connected.degree[node] == 1]
        if not endpoints:
            endpoints = list(comp)  # fallback: gebruik alle knopen
        component_endpoints.append(endpoints)

    # Maak lijst van alle mogelijke verbindingen tussen componenten
    candidate_edges = []
    for i, j in combinations(range(len(components)), 2):
        for u in component_endpoints[i]:
            for v in component_endpoints[j]:
                d = dist(u, v)
                candidate_edges.append((d, u, v))

    # Sorteer op afstand
    candidate_edges.sort()

    # Gebruik Union-Find om componenten te verbinden
    parent = {node: node for comp in components for node in comp}

    def find(u):
        while parent[u] != u:
            parent[u] = parent[parent[u]]
            u = parent[u]
        return u

    def union(u, v):
        parent[find(u)] = find(v)

    # Voeg kortste verbindingen toe totdat alles verbonden is
    for d, u, v in candidate_edges:
        if find(u) != find(v):
            G_connected.add_edge(u, v, weight=d)
            union(u, v)

    return G_connected


# def connect_disconnected_networks(multilinestring):
#     """
#     Takes a MultiLineString consisting of multiple disconnected networks and returns
#     a new MultiLineString where the minimal set of connecting lines is added to form
#     one connected network.
#     """
#     # Merge lines and extract connected components
#     merged = linemerge(multilinestring)
#     if isinstance(merged, LineString):
#         return merged  # Already connected
#
#     # Build graph from lines
#     G = nx.Graph()
#     for line in merged:
#         coords = list(line.coords)
#         G.add_edge(coords[0], coords[-1], geometry=line)
#
#     # Find connected components
#     components = list(nx.connected_components(G))
#     if len(components) <= 1:
#         return merged  # Already connected
#
#     # Get representative points (endpoints) from each component
#     component_endpoints = []
#     for comp in components:
#         endpoints = []
#         for node in comp:
#             if G.degree[node] == 1:
#                 endpoints.append(Point(node))
#         component_endpoints.append(endpoints)
#
#     # Build a list of candidate connections between components
#     connection_lines = []
#     connected = set()
#     remaining = set(range(len(component_endpoints)))
#
#     # Use a greedy approach to connect components
#     while len(remaining) > 1:
#         min_dist = math.inf
#         best_pair = None
#         best_line = None
#         for i in remaining:
#             for j in remaining:
#                 if i >= j:
#                     continue
#                 for p1 in component_endpoints[i]:
#                     for p2 in component_endpoints[j]:
#                         dist = p1.distance(p2)
#                         if dist < min_dist:
#                             min_dist = dist
#                             best_pair = (i, j)
#                             best_line = LineString([p1, p2])
#         if best_pair:
#             connection_lines.append(best_line)
#             remaining.remove(best_pair[1])
#             # Merge endpoints of j into i
#             component_endpoints[best_pair[0]].extend(component_endpoints[best_pair[1]])
#
#     # Combine original and connection lines
#     all_lines = list(merged) + connection_lines
#     return MultiLineString(all_lines)


def prepare_network(segments):
    network = line_merge(safe_unary_union(segments))
    network = snap(network, network, tolerance=0.01)
    # TODO check if this is sufficient; possibly we also need to inject vertex on edge where line almost touches an edge
    network = fill_gaps_in_multilinestring(
        network, 0.1
    )  # also needed to fill 'gaps' to connect reference objects fe points
    return line_merge(safe_unary_union(network))


def insert_vertex(multilinestring, point):
    new_lines = []
    if isinstance(multilinestring, LineString):
        multilinestring = to_multi(multilinestring, None)

    for line in multilinestring.geoms:
        # Find the nearest point on the line to the given point
        nearest_point = line.interpolate(line.project(point))

        # Split the line at the nearest point
        coords = list(line.coords)
        min_dist = float("inf")
        insert_index = None

        for i in range(len(coords) - 1):
            segment = LineString([coords[i], coords[i + 1]])
            dist = segment.distance(nearest_point)
            if dist < min_dist:
                min_dist = dist
                insert_index = i + 1

        # Insert the nearest point into the coordinates
        new_coords = (
            coords[:insert_index]
            + [(nearest_point.x, nearest_point.y)]
            + coords[insert_index:]
        )
        new_lines.append(LineString(new_coords))

    return MultiLineString(new_lines)


def total_vertex_distance(
    geom1: BaseGeometry, geom2: BaseGeometry, bidirectional=True
) -> float:
    """Compute the total vertex-based distance between two geometries."""
    vertices1 = get_coords_from_geometry(geom1)
    total_distance_1 = 0.0
    len_vertices_1 = len(vertices1)
    # Sum of distances from vertices1 to geom2
    for pt in vertices1:
        total_distance_1 += Point(pt).distance(geom2)
    total_distance = total_distance_1
    len_vertices = len_vertices_1

    if bidirectional:
        # Sum of distances from vertices2 to geom1
        total_distance_2 = 0.0
        vertices2 = get_coords_from_geometry(geom2)
        len_vertices_2 = len(vertices2)
        for pt in vertices2:
            total_distance_2 += Point(pt).distance(geom1)
        total_distance = total_distance + total_distance_2
        len_vertices = len_vertices + len_vertices_2

    return total_distance / len_vertices


import networkx as nx
from shapely.geometry import LineString, Point
import math


def euclidean_distance(p1, p2):
    return math.hypot(p1[0] - p2[0], p1[1] - p2[1])


def hausdorff_distance(ls1, ls2):
    return ls1.hausdorff_distance(ls2)


def path_to_linestring(path):
    return LineString(path)


def find_cycle_paths(cycle):
    """Split a cycle into two alternative paths between the furthest nodes"""
    max_dist = 0
    start, end = cycle[0], cycle[1]
    for i in range(len(cycle)):
        for j in range(i + 1, len(cycle)):
            d = euclidean_distance(cycle[i], cycle[j])
            if d > max_dist:
                max_dist = d
                start, end = cycle[i], cycle[j]

    idx_start = cycle.index(start)
    idx_end = cycle.index(end)

    if idx_start < idx_end:
        path1 = cycle[idx_start : idx_end + 1]
        path2 = cycle[idx_end:] + cycle[: idx_start + 1]
    else:
        path1 = cycle[idx_end : idx_start + 1]
        path2 = cycle[idx_start:] + cycle[: idx_end + 1]

    return path1, path2


def remove_path_from_graph(G, path):
    for i in range(len(path) - 1):
        if G.has_edge(path[i], path[i + 1]):
            G.remove_edge(path[i], path[i + 1])


def simplify_graph_by_best_cycle_path(G, input_line):
    """
    Vereenvoudigt een undirected graph door in elke cycle enkel het pad te behouden
    dat het dichtst ligt bij een opgegeven LineString.

    Parameters:
    - G: NetworkX Graph met nodes die een 'pos' attribuut hebben (x, y)
    - input_line: Shapely LineString waartegen afstand gemeten wordt

    Returns:
    - G_simplified: een nieuwe NetworkX Graph met vereenvoudigde edges
    """
    # Zet om naar directed graph voor cycle detectie
    DG = nx.DiGraph(G)
    cycles = list(nx.simple_cycles(DG))
    edges_to_keep = set()

    for cycle in cycles:
        # Genereer beide richtingen van het pad in de cycle
        paths = [
            [(cycle[i], cycle[(i + 1) % len(cycle)]) for i in range(len(cycle))],
            [(cycle[i], cycle[i - 1]) for i in range(len(cycle))],
        ]

        path_distances = []

        for path in paths:
            total_distance = 0
            for u, v in path:
                p1 = Point(u)
                p2 = Point(v)
                edge_line = LineString([p1, p2])
                total_distance += edge_line.distance(input_line)
            avg_distance = total_distance / len(path)
            path_distances.append((avg_distance, path))

        # Kies het pad met de kleinste gemiddelde afstand
        best_path = min(path_distances, key=lambda x: x[0])[1]
        edges_to_keep.update(best_path)

    # Bouw de vereenvoudigde graaf
    G_simplified = nx.Graph()
    for node, data in G.nodes(data=True):
        G_simplified.add_node(node, **data)

    for u, v in G.edges():
        if (u, v) in edges_to_keep or (v, u) in edges_to_keep:
            G_simplified.add_edge(u, v)

    return G_simplified


def find_best_path_in_network(
    geom_to_process, nw_multilinestring, snap_strategy, relevant_distance
):
    """
    Detrlmine the best path between 2 points in the network using the Hausdorf-distance
    Parameters:
    - geom_to_process: shapely.geometry.MultiLineString -  geometry (line) with startpoint-endpoint and to determine hausdorf-distance
    - nw_multilinestring: shapely.geometry.MultiLineString - MultiLineString with all parts of a network

    Returns:
    - shapely.geometry.LineString van het langste pad tussen de twee punten
    """
    # if isinstance(nw_multilinestring, LineString):
    #     return nw_multilinestring

    start_point = Point(geom_to_process.coords[0])
    end_point = Point(geom_to_process.coords[-1])
    nw_multilinestring = insert_vertex(nw_multilinestring, start_point)
    nw_multilinestring = insert_vertex(nw_multilinestring, end_point)

    # Create graph
    G = graph_from_multilinestring(nw_multilinestring)
    # remove cycles #todo test if this improves
    # G,removed_edges = simplify_graph_by_best_cycle_path(G,geom_to_process)

    start_node = nearest_node(start_point, G.nodes)
    end_node = nearest_node(end_point, G.nodes)

    if start_node == end_node:
        return find_best_circle_path(G.to_directed(), geom_to_process)

    if not snap_strategy is None and snap_strategy != SnapStrategy.NO_PREFERENCE:
        # when SnapStrategy = PREFER_VERTICES or ONLY_VERTICES
        start_node = get_vertex_node(G, relevant_distance, start_node, start_point)
        end_node = get_vertex_node(G, relevant_distance, end_node, end_point)

    # Search all simple paths (limited to 100, because cyclic paths can result in a lot of simple paths
    all_paths = islice(nx.all_simple_paths(G, source=start_node, target=end_node), 1000)

    # Determine the network-path that fits the best to the original inputgeometry
    min_dist = inf
    best_line = None
    # i=0
    for path in all_paths:
        # i=i+1
        # print(str(i))
        line = LineString(path)
        dist = total_vertex_distance(line, geom_to_process)
        if dist < min_dist:
            min_dist = dist
            best_line = line
    return best_line


def get_vertex_node(G, relevant_distance, input_node, point):
    min_dist = inf
    for node in G.neighbors(input_node):
        dist = point.distance(Point(node)) * BUFFER_MULTIPLICATION_FACTOR
        # TODO is there a better way to check which nodes to use? If the distance is almost the same as the relevant distance, it could be that the vertex is constructed vertex by reference_intersection
        if dist < relevant_distance and dist < min_dist:
            min_dist = dist
            input_node = node
    return input_node


def remove_pseudonodes(G):
    # TODO; research if removing pseudonodes can improve performance.
    # Maybe better to control when making the initial graph?
    G = G.copy()
    for node in list(G.nodes):
        if G.degree[node] == 2:
            neighbors = list(G.neighbors(node))
            if len(neighbors) == 2:
                # Voeg een nieuwe edge toe tussen de buren
                if not G.has_edge(neighbors[0], neighbors[1]):
                    G.add_edge(neighbors[0], neighbors[1])
                G.remove_node(node)
    return G


# def scale_segments(segments, factor=1.01):
#     scaled_segments = []
#     for segment in segments:
#         scaled_segment = scale(segment, xfact=factor, yfact=factor, origin="center")
#         scaled_segments.append(scaled_segment)
#     return scaled_segments


# def extend_line_in_graph(existing_line, graph):
#     # TODO; check if we are going to use this
#
#     start_node = nearest_node(Point(existing_line.coords[0]), graph.nodes)
#     end_node = nearest_node(Point(existing_line.coords[-1]), graph.nodes)
#
#     # Zoek de bijbehorende node in de graph
#     # (je moet hier mogelijk een mapping hebben tussen coördinaten en node IDs)
#
#     start_geom = get_next_geometry(graph, start_node)
#     end_geom = get_next_geometry(graph, end_node)
#
#     # Voeg toe aan bestaande lijn
#     return safe_unary_union([existing_line, start_geom, end_geom])
#
#
# def get_next_geometry(graph, node):
#     # TODO; check if we are going to use this
#
#     # Stel dat je de node ID hebt:
#     neighbors = list(graph.neighbors(node))
#     # Kies een volgende node (bijv. de eerste)
#     next_node = neighbors[0]
#     # Haal de edge geometrie op
#     edge_data = graph.get_edge_data(node, next_node)
#     return edge_data[0]["geometry"]


def _get_snapped_point(point, ref_line, ref_coords, snap_strategy, tolerance):
    p1, p2 = None, None
    p1_vertices, p2_vertices = None, None
    distance_ref_line = inf
    if not ref_line.is_empty:
        p1, p2 = nearest_points(point, ref_line)
        distance_ref_line = p2.distance(point)
    if len(ref_coords) != 0:
        p1_vertices, p2_vertices = nearest_points(point, MultiPoint(ref_coords))
        distance_ref_coords = p2_vertices.distance(point)
        if distance_ref_coords < distance_ref_line:
            p1, p2 = p1_vertices, p2_vertices
    return _snapped_point_by_snapstrategy(
        point, p2, p2_vertices, snap_strategy, tolerance
    )


def add_point_as_node_on_closest_edge(G, point):
    """
    Voeg een Shapely Point toe als node op de dichtstbijzijnde edge in een NetworkX-graaf.
    De originele edge wordt gesplitst in twee nieuwe edges, met aangepaste lengtes.

    Parameters:
    - G: NetworkX-graaf met nodes als (x, y)-tuples
    - point: Shapely Point

    Returns:
    - G: aangepaste graaf met nieuwe node en gesplitste edges
    """
    min_distance = float("inf")
    closest_edge = None
    projected_point = None

    # Zoek de dichtstbijzijnde edge
    for u, v, data in G.edges(data=True):
        line = LineString([u, v])
        proj = line.interpolate(line.project(point))
        distance = point.distance(proj)
        if distance < min_distance:
            min_distance = distance
            closest_edge = (u, v)
            projected_point = proj

    if closest_edge is None or projected_point is None:
        raise ValueError("Geen geschikte edge gevonden.")

    u, v = closest_edge
    original_length = G[u][v].get("length", dist(u, v))
    G.remove_edge(u, v)

    new_node = (projected_point.x, projected_point.y)
    G.add_node(new_node)

    # Bereken nieuwe lengtes
    length1 = dist(u, new_node)
    length2 = dist(new_node, v)

    G.add_edge(u, new_node, length=length1)
    G.add_edge(new_node, v, length=length2)

    return G


# def multilinestring_to_graph(mls):
#     """
#     Converts a Shapely MultiLineString into a NetworkX graph.
#     Nodes are only endpoints and intersection points.
#     Edges represent segments between these nodes and include length as an attribute.
#     """
#     # Step 1: Merge all lines into a single geometry and find intersections
#     merged = linemerge(mls)
#     if isinstance(merged, LineString):
#         lines = [merged]
#     else:
#         lines = list(merged)
#
#     # Step 2: Collect all endpoints
#     endpoints = set()
#     for line in lines:
#         endpoints.add(Point(line.coords[0]))
#         endpoints.add(Point(line.coords[-1]))
#
#     # Step 3: Find all intersection points
#     intersections = set()
#     for i, line1 in enumerate(lines):
#         for j, line2 in enumerate(lines):
#             if i < j:
#                 inter = line1.intersection(line2)
#                 if "Point" in inter.geom_type:
#                     intersections.add(inter)
#                 elif inter.geom_type == "MultiPoint":
#                     intersections.update(inter.geoms)
#
#     # Combine endpoints and intersections
#     split_points = list(endpoints.union(intersections))
#
#     # Step 4: Split lines at split_points
#     split_lines = []
#     for line in lines:
#         for pt in split_points:
#             if not line.contains(pt):
#                 continue
#             line = snap(line, pt, 1e-8)
#         result = split(line, unary_union(split_points))
#         split_lines.extend(result.geoms)
#
#     # Step 5: Build graph
#     G = nx.Graph()
#     for segment in split_lines:
#         coords = list(segment.coords)
#         if len(coords) < 2:
#             continue
#         p1 = tuple(coords[0])
#         p2 = tuple(coords[-1])
#         length = segment.length
#         G.add_edge(p1, p2, length=length)
#
#     return G


def extract_points_lines_from_geometry(geometry):
    """
    Recursief haalt deze functie alle Points en LineStrings uit een Shapely-geometrie.
    Polygonen worden omgezet naar hun exterior en interior ringen als LineStrings.
    GeometryCollections worden volledig doorlopen.
    """
    geometries = []

    if geometry.is_empty:
        return GeometryCollection()

    if isinstance(geometry, Polygon):
        geometries.append(LineString(geometry.exterior.coords))
        for interior in geometry.interiors:
            geometries.append(LineString(interior.coords))

    elif isinstance(geometry, LineString):
        geometries.append(geometry)

    elif isinstance(
        geometry, (MultiPoint, MultiPolygon, MultiLineString, GeometryCollection)
    ):
        for geom in geometry.geoms:
            geometries.extend(extract_points_lines_from_geometry(geom).geoms)

    elif isinstance(geometry, (Point)):
        geometries.append(geometry)

    return GeometryCollection(geometries)
