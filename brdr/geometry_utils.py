import logging
from itertools import combinations
from math import pi, inf

import networkx as nx
import numpy as np
from shapely import GEOSException, equals, shortest_line
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
from shapely import symmetric_difference
from shapely import to_wkt
from shapely import unary_union
from shapely import union
from shapely.affinity import scale
from shapely.geometry.base import BaseGeometry
from shapely.lib import line_merge
from shapely.ops import substring, nearest_points
from shapely.prepared import prep

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
    # TODO - cleanup; this gives errors with Geometrycollection-empty
    result = remove_shortest_and_merge(
        result, relevant_length=0.01
    )  # Solves small (0.01) overshoots in multilinestrings #TODO; maybe this has to be cleaned in another place?
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


from shapely.ops import snap


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


# Stap 2: vind de dichtstbijzijnde knopen bij start- en eindpunt
def nearest_node(point, nodes):
    return min(nodes, key=lambda n: Point(n).distance(point))


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


def find_longest_path_between_points(multilinestring, start_point, end_point):
    """
    Bepaal het langste pad tussen twee punten in een MultiLineString-geometrie.

    Parameters:
    - multilinestring: shapely.geometry.MultiLineString
    - start_point: shapely.geometry.Point
    - end_point: shapely.geometry.Point

    Returns:
    - shapely.geometry.LineString van het langste pad tussen de twee punten
    """
    if isinstance(multilinestring, LineString):
        return multilinestring
    # Stap 1: maak een graaf van de lijnen
    G = nx.Graph()
    for line in multilinestring.geoms:
        coords = list(line.coords)
        for i in range(len(coords) - 1):
            p1 = tuple(coords[i])
            p2 = tuple(coords[i + 1])
            segment = LineString([p1, p2])
            G.add_edge(p1, p2, weight=segment.length, geometry=segment)

    start_node = nearest_node(start_point, G.nodes)
    end_node = nearest_node(end_point, G.nodes)
    if start_node == end_node:
        return find_circle_path(G.to_directed())

    # Stap 3: zoek alle eenvoudige paden tussen start en eind
    all_paths = nx.all_simple_paths(G, source=start_node, target=end_node)

    # Stap 4: bepaal het pad met de grootste som van gewichten
    max_length = 0
    longest_path = []
    for path in all_paths:
        length = sum(G[path[i]][path[i + 1]]["weight"] for i in range(len(path) - 1))
        if length > max_length:
            max_length = length
            longest_path = path

    # Stap 5: zet het pad om naar een LineString
    if longest_path:
        return LineString(longest_path)
    else:
        return None


def longest_path_from_multilinestring(mls):
    if isinstance(mls, LineString):
        return mls
    # Create an undirected graph
    G = nx.Graph()

    # Add edges to the graph from the MultiLineString
    for line in mls.geoms:
        coords = list(line.coords)
        for i in range(len(coords) - 1):
            start = coords[i]
            end = coords[i + 1]
            length = LineString([start, end]).length
            G.add_edge(start, end, weight=length)

    # Find the longest simple path by checking all simple paths

    longest_path = []
    max_length = 0

    for source in G.nodes:
        for target in G.nodes:
            if source != target:
                try:
                    for path in nx.all_simple_paths(G, source=source, target=target):
                        path_length = sum(
                            LineString([path[i], path[i + 1]]).length
                            for i in range(len(path) - 1)
                        )
                        if path_length > max_length:
                            max_length = path_length
                            longest_path = path
                except nx.NetworkXNoPath:
                    continue

    # Convert the longest path to a LineString
    if longest_path:
        return LineString(longest_path)
    else:
        return None


# def longest_connected_linestring(multilinestring):
#     # Merge connected lines into LineStrings
#     merged = line_merge(unary_union(multilinestring))
#
#     # Ensure result is iterable (could be a single LineString or MultiLineString)
#     if isinstance(merged, LineString):
#         return merged
#     elif isinstance(merged, MultiLineString):
#         # Return the longest LineString
#         return max(merged.geoms, key=lambda line: line.length)
#     else:
#         return None


def scale_segments(segments, factor=1.01):
    scaled_segments = []
    for segment in segments:
        scaled_segment = scale(segment, xfact=factor, yfact=factor, origin="center")
        scaled_segments.append(scaled_segment)
    return scaled_segments


def remove_shortest_and_merge(multilinestring, relevant_length=float("inf")):
    """
    Tries to merge a multilinestring to a linestring, and recursively tries to remove the shortest part of a multilinestring until it can be merged to a linestring,
    :param multilinestring: MultiLineString to merge
    :param relevant_length: Value when the length of a part of multilinestring is long enough to be relevant and multilinestring is kept. Default is infinity so it will recursively end into a single LineString.    :return:
    """
    if multilinestring.is_empty:
        return multilinestring
    # Check if the merged result is a LineString
    multilinestring = line_merge(safe_unary_union(multilinestring))
    if isinstance(multilinestring, LineString):
        return multilinestring

    # Find the shortest LineString in the MultiLineString
    lines = list(multilinestring.geoms)
    distance_list = [line.length for line in lines]
    shortest_distance = min(distance_list)  # find the line closest to the point
    shortest_line = lines[distance_list.index(shortest_distance)]
    if shortest_distance >= relevant_length:
        return multilinestring
    # Create a new MultiLineString without the shortest LineString
    remaining_lines = [line for line in lines if line != shortest_line]
    new_multilinestring = MultiLineString(remaining_lines)
    # Recursively call the function with the new MultiLineString
    return remove_shortest_and_merge(new_multilinestring)


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


from shapely.geometry import (
    Point,
    LineString,
    Polygon,
    MultiPoint,
    MultiLineString,
    MultiPolygon,
    GeometryCollection,
)


def extract_points_lines_from_geometry(geometry):
    """
    Recursief haalt deze functie alle LineStrings uit een Shapely-geometrie.
    Polygonen worden omgezet naar hun exterior en interior ringen als LineStrings.
    GeometryCollections worden volledig doorlopen.
    Andere geometrieën zoals punten worden genegeerd.
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
