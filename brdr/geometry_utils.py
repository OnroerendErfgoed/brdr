import logging
import math
from itertools import combinations
from math import inf
from math import pi
from typing import Union, List, Tuple

import networkx as nx
import numpy as np
import pyproj
from networkx import Graph
from shapely import GEOSException, get_coordinates, LinearRing
from shapely import STRtree
from shapely import buffer
from shapely import difference
from shapely import equals
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
from shapely.errors import GeometryTypeError
from shapely.geometry import (
    Point,
    LineString,
    Polygon,
    MultiPoint,
    MultiLineString,
    MultiPolygon,
    GeometryCollection,
    base,
)
from shapely.geometry.base import BaseGeometry
from shapely.lib import line_merge
from shapely.ops import nearest_points, split
from shapely.ops import substring
from shapely.prepared import prep

from brdr.enums import SnapStrategy

ShapelyGeometry = Union[
    Point,
    LineString,
    Polygon,
    MultiPoint,
    MultiLineString,
    MultiPolygon,
    GeometryCollection,
]

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
        >>> logging.debug(result)
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
        >>> logging.debug(result)
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
        >>> logging.debug(result)
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
        >>> logging.debug(result)
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
        unioned_geom = make_valid(unary_union(geometries))
    except:
        geometries = [make_valid(g) for g in geometries]
        unioned_geom = make_valid(unary_union(geometries))

    return unioned_geom


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
    ref_vertices_objs = [Point(c) for c in ref_coords]
    tree = STRtree(ref_vertices_objs) if ref_vertices_objs else None

    if len(ref_coords) == 0:
        snap_strategy = SnapStrategy.NO_PREFERENCE

    points = []
    for geom in geometry.geoms:
        coords = list(geom.coords)
        coordinates = []
        for idx, coord in enumerate(coords):
            p = Point(coords[idx])
            p_snapped, bool_snapped, ref_vertices = _get_snapped_point(
                p, ref_border, tree, ref_vertices_objs, snap_strategy, tolerance
            )
            coordinates.append(p_snapped.coords[0])

        # convert coordinates back to a point
        point = make_valid(Point(coordinates))
        points.append(point)
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
        lines.append(linestring)
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
        polygons.append(polygon)
    return buffer_neg_pos(safe_unary_union(polygons), correction_distance)


def _get_ref_objects(reference):
    # reference = safe_intersection(
    #     reference, buffer_pos(geometry, tolerance * BUFFER_MULTIPLICATION_FACTOR)
    # )
    ref_coords = list(get_coordinates(reference))
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
    coords: List[Tuple[float, float]],
    ref_border: LineString,
    ref_borders: List[LineString],
    ref_coords: List[Tuple[float, float]],
    snap_strategy: str,
    tolerance: float,
) -> List[Tuple[float, float]]:
    """
    Snaps a sequence of coordinates to reference geometries using optimized spatial lookups.

    This function iterates through vertices of a polyline and snaps them to a
    reference border or specific vertices. It optimizes performance by using an
    STRtree for nearest-neighbor lookups and pre-calculating buffers to avoid
    redundant geometric operations within the loop.

    Parameters
    ----------
    coords : list of tuple of float
        The input coordinates (vertices) to be snapped.
    ref_border : LineString
        The primary reference line used for snapping.
    ref_borders : list of LineString
        A collection of reference borders used for sub-linestring interpolation.
    ref_coords : list of tuple of float
        A list of (x, y) coordinates representing reference vertices.
    snap_strategy : str
        The snapping strategy to employ (e.g., 'PREFER_VERTICES').
    tolerance : float
        The maximum distance allowed for a coordinate to be snapped.

    Returns
    -------
    list of tuple of float
        A list of snapped coordinates, including interpolated vertices where
        snapping to a reference path was successful.
    """
    if not coords:
        return []

    coordinates = []

    # 1. Pre-process references (Outside the loop!)
    # Pre-buffering avoids thousands of expensive re-calculations inside the loop
    ref_buffered_main = buffer_pos(ref_border, tolerance)
    ref_coords_multipoint = MultiPoint(ref_coords)
    ref_coords_buffered = buffer_pos(ref_coords_multipoint, tolerance)

    # Initialize Spatial Index (STRtree) for log(N) vertex lookups
    ref_vertices_objs = [Point(c) for c in ref_coords]
    tree = STRtree(ref_vertices_objs) if ref_vertices_objs else None

    # 2. Cache the snapped start point for the first iteration
    p_start = Point(coords[0])
    p_start_snapped, bool_start_snapped, ref_vertices_start = _get_snapped_point(
        p_start, ref_border, tree, ref_vertices_objs, snap_strategy, tolerance
    )

    for idx in range(1, len(coords)):
        p_end = Point(coords[idx])

        # Only snap the end point (the start point is carried over from the previous iteration)
        p_end_snapped, bool_end_snapped, ref_vertices_end = _get_snapped_point(
            p_end, ref_border, tree, ref_vertices_objs, snap_strategy, tolerance
        )

        coordinates.append(p_start_snapped.coords[0])

        if not bool_start_snapped and not bool_end_snapped:
            # No snapping required for either point, proceed normally
            coordinates.append(p_end_snapped.coords[0])

        elif bool_start_snapped and bool_end_snapped:
            # Both points snapped: retrieve the path along the reference borders
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

        else:
            # Mixed scenario: one point is snapped, the other is not
            line = LineString([p_start, p_end])

            # Select the appropriate pre-buffered geometry
            curr_ref_buffered = (
                ref_coords_buffered
                if (ref_vertices_start or ref_vertices_end)
                else ref_buffered_main
            )
            intersected_line = safe_intersection(line, curr_ref_buffered)

            if (
                not intersected_line.is_empty
                and intersected_line.geom_type == "LineString"
            ):
                boundary = intersected_line.boundary.geoms
                first, last = boundary[0], boundary[-1]

                # Efficiently determine the midpoint for snapping
                p_mid = last if first == p_start else first if last == p_end else None

                if p_mid:
                    p_mid_snapped, _, _ = _get_snapped_point(
                        p_mid,
                        ref_border,
                        tree,
                        ref_vertices_objs,
                        snap_strategy,
                        tolerance,
                    )

                    if first == p_start:
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
                    else:
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

        # 3. Prepare for next iteration: current 'end' becomes next 'start'
        p_start = p_end
        p_start_snapped = p_end_snapped
        bool_start_snapped = bool_end_snapped
        ref_vertices_start = ref_vertices_end

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
    for line in list(to_multi(line_substring).geoms):
        for p in line.coords:
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
    ##259 Research around optimization of partitioning
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
    return str(geometry.bounds).strip("()").replace(" ", "")


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


# def get_geoms_from_geometry(geometry):
#     geoms = set()
#     if geometry is None or geometry.is_empty:
#         return geoms
#     elif isinstance(geometry, (Point, LineString, Polygon)):
#         geoms.update([geometry])
#     elif isinstance(geometry, (MultiPoint, MultiLineString, MultiPolygon)):
#         geoms.update(geometry.geoms)
#     elif isinstance(geometry, GeometryCollection):
#         for geom in geometry.geoms:
#             geoms.update(get_geoms_from_geometry(geom))
#     return geoms


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
    find closest node to point
    :param point:
    :param nodes:
    :return:
    """
    return min(nodes, key=lambda n: Point(n).distance(point))


# def find_best_circle_path2(G, geom_to_process):
#     """
#             #TODO - research if this is a better alternative for 'find_cycle_path'?
#     """
#     line_geoms = []
#     for edge in G.edges(data=True):
#         line_geoms.append(edge[2]["geometry"])
#     all_faces = polygonize(line_geoms)
#     all_polygons = [geom for geom in all_faces.geoms if isinstance(geom, Polygon)]
#     geom_to_process_polygon = Polygon(geom_to_process)
#     candidates = [face for face in all_polygons if face.intersects(geom_to_process_polygon)]
#
#     best_iou = 0
#
#     must_haves = []
#
#     remaining_data = []  # Lijst van tuples: (polygoon, overlap_area)
#
#     # Eén pass om te splitsen
#     for candidate in candidates:
#         c_area = candidate.area
#         i_area = candidate.intersection(geom_to_process_polygon).area
#
#         if i_area >= c_area * 0.99:
#             must_haves.append(candidate)
#         elif i_area > 0:  # Alleen meenemen als er überhaupt overlap is
#             # We slaan de i_area even op in het object om niet opnieuw te hoeven rekenen
#             remaining_data.append((candidate, i_area))
#
#     # Directe start met de bulk
#     current_selection = must_haves
#     combined_union = unary_union(must_haves) if must_haves else None
#
#     # 2. Sorteer de twijfelgevallen op de overlap-waarde (index 1 van de tuple)
#     remaining_data.sort(key=lambda x: x[1], reverse=True)
#
#     # Initialiseer best_iou
#     if combined_union:
#         intersection = combined_union.intersection(geom_to_process_polygon).area
#         union = combined_union.union(geom_to_process_polygon).area
#         best_iou = intersection / union if union > 0 else 0
#     else:
#         best_iou = 0
#
#     # 3. De Greedy loop over de uitgepakte tuples
#     for candidate, i_area in remaining_data:
#         test_set = current_selection + [candidate]
#         test_combined = unary_union(test_set)
#
#         test_intersection = test_combined.intersection(geom_to_process_polygon).area
#         test_union = test_combined.union(geom_to_process_polygon).area
#         iou = test_intersection / test_union if test_union > 0 else 0
#
#         if iou > best_iou:
#             best_iou = iou
#             current_selection.append(candidate)
#         else:
#             continue
#     try:
#         final_cycle = unary_union(current_selection)
#         if final_cycle is None or final_cycle.is_empty:
#             return None
#
#         final_cycle=to_multi(final_cycle)
#         final_cycle = [geom for geom in final_cycle.geoms if isinstance(geom, Polygon)]
#         final_cycle = unary_union(final_cycle)
#         outer_line = final_cycle.exterior
#     except:
#         outer_line = None
#         #raise Exception
#     return outer_line

def find_best_circle_path(graph, geom_to_process):
    """
    Find the best cycle in the graph that is closest to the original geometry to process
    """
    cycles_generator = nx.cycle_basis(graph)
    min_dist = inf
    best_cycle_line = None
    max_amount = 1000

    for i, cycle in enumerate(cycles_generator):
        if len(cycle) <= 2:
            continue
        if i > max_amount:  # safetyleak
            logging.warning(
                f"max cycles tested while searching for geometry: {max_amount}"
            )
            break
        # logging.debug (i)
        # logging.debug (LineString(cycle).wkt)
        cycle_coords = cycle + [cycle[0]]

        cycle_line = LineString(cycle_coords)
        vertex_distance = total_vertex_distance(cycle_line, geom_to_process)
        if vertex_distance < min_dist:
            min_dist = vertex_distance
            best_cycle_line = cycle_line
    return best_cycle_line


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
    graph = Graph()
    _multilinestring_to_edges(graph, multilinestring, node_tag="",edge_tag="",pseudo_coords=[],pseudonode_tag="",)

    # Find all simple paths and keep the longest one
    longest_path = []
    max_length = 0

    for source in graph.nodes:
        for target in graph.nodes:
            if source != target:
                try:
                    path = nx.shortest_path(
                        graph, source=source, target=target, weight="weight"
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


def euclidean_distance(p1, p2):
    return math.hypot(p1[0] - p2[0], p1[1] - p2[1])


def total_vertex_distance(
    geom1: BaseGeometry, geom2: BaseGeometry, bidirectional=True
) -> float:
    """Compute the total vertex-based distance between two geometries."""
    vertices1 = get_coordinates(geom1)
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
        vertices2 = get_coordinates(geom2)
        len_vertices_2 = len(vertices2)
        for pt in vertices2:
            total_distance_2 += Point(pt).distance(geom1)
        total_distance = total_distance + total_distance_2
        len_vertices = len_vertices + len_vertices_2

    return total_distance / len_vertices


def remove_composite_edges(G: nx.Graph) -> List[Tuple]:
    """
    Identifies and removes edges that are exact compositions of other edges.

    An edge (u, v) is considered composite if there exists an alternative
    path between u and v whose total weight is exactly equal to the weight
    of the direct edge.

    Parameters
    ----------
    G : nx.Graph
        The undirected graph to be processed. This object is modified in-place.

    Returns
    -------
    removed_edges : list of tuple
        A list of edges (u, v) that were removed from the graph.
    """
    removed_edges = []

    # Iterate over a static list of edges to avoid 'dictionary changed size' errors
    for u, v in list(G.edges()):
        # Get the weight of the current direct edge
        direct_weight = G[u][v].get("weight", 1.0)

        # Temporarily remove the edge to find alternative paths
        G.remove_edge(u, v)

        try:
            # Find the shortest path using the remaining edges
            alt_path_len = nx.shortest_path_length(G, u, v, weight="weight")

            # If an alternative path exists with the EXACT same weight,
            # the direct edge is a composition and remains removed.
            if np.isclose(alt_path_len, direct_weight):
                removed_edges.append((u, v))
            else:
                # Path is longer or shorter (shortcut), so we restore the edge
                G.add_edge(u, v, weight=direct_weight)

        except nx.NetworkXNoPath:
            # No alternative route exists, restore the direct edge
            G.add_edge(u, v, weight=direct_weight)

    return removed_edges


def find_best_path_in_network(geom_to_process, graph):
    """
    Determine the best path between 2 points in the network using the Hausdorf-distance
    Parameters:
    - geom_to_process: shapely.geometry.MultiLineString -  geometry (line) with startpoint-endpoint and to determine hausdorf-distance
    - nw_multilinestring: shapely.geometry.MultiLineString - MultiLineString with all parts of a network

    Returns:
    - shapely.geometry.LineString van het langste pad tussen de twee punten
    """

    start_point = Point(geom_to_process.coords[0])
    end_point = Point(geom_to_process.coords[-1])

    start_node = nearest_node(start_point, graph.nodes)
    end_node = nearest_node(end_point, graph.nodes)

    if start_node == end_node:
        return find_best_circle_path(graph, geom_to_process)

    # Search all simple paths (limited, because cyclic paths can result in a lot of simple paths
    max_amount = 1000
    path_found = nx.has_path(graph, start_node, end_node)
    logging.debug("Path detected? - " + str(path_found))
    if not path_found:
        graph = finalize_network_connectivity(
            graph,
            50,
            edge_tags=["theme_lines", "ref_lines", "interconnect", "gap_closure"],
        )

    all_paths_generator = nx.all_simple_paths(
        graph, source=start_node, target=end_node, cutoff=max_amount
    )

    # Determine the network-path that fits the best to the original inputgeometry
    min_dist = inf
    best_line = None
    for i, path in enumerate(all_paths_generator):
        if i > max_amount:  # safetyleak
            logging.warning(
                f"max paths tested while searching for geometry: {max_amount}"
            )
            break
        try:  # added try/except because 'path' sometimes exists out of 1 point, resulting in LineString-error
            line = LineString(path)
            dist = total_vertex_distance(line, geom_to_process)
            if dist < min_dist:
                min_dist = dist
                best_line = line
        except:
            pass
    return best_line

def _get_snapped_point(
    point: Point,
    ref_line: LineString,
    ref_coords_tree: STRtree,
    ref_vertices_objs: List[Point],
    snap_strategy: str,
    tolerance: float,
) -> Tuple[Point, bool, bool]:
    """
    Determines the nearest point to the input from a reference line or spatial tree.

    This function optimizes the snapping process by using a pre-computed STRtree
    for vertex lookups and compares distances between the nearest point on a
    line versus the nearest vertex to determine the best snapping candidate.

    Parameters
    ----------
    point : Point
        The source point to be snapped.
    ref_line : LineString
        The reference line geometry to snap onto.
    ref_coords_tree : STRtree
        A spatial index tree containing the reference vertices for fast lookup.
    ref_vertices_objs : list of Point
        The original list of Point objects used to build the STRtree,
        required for index-based lookups in certain Shapely versions.
    snap_strategy : str
        The snapping logic to apply (e.g., 'PREFER_VERTICES', 'ONLY_VERTICES').
    tolerance : float
        The maximum distance allowed for a successful snap.

    Returns
    -------
    p_snapped : Point
        The resulting point (either the original, or the snapped location).
    snapped : bool
        True if the point was moved/snapped within the tolerance.
    is_vertex : bool
        True if the point was snapped specifically to a vertex.
    """
    p_nearest_line = None
    p_nearest_vertex = None
    dist_line = float("inf")
    dist_vertex = float("inf")

    # 1. Check proximity to the reference line
    if not ref_line.is_empty:
        # p_on_line is the projection of the point onto the line
        p_on_line, _ = nearest_points(ref_line, point)
        dist_line = p_on_line.distance(point)
        p_nearest_line = p_on_line

    # 2. Check proximity to vertices using the STRtree index
    if ref_coords_tree is not None:
        result = ref_coords_tree.nearest(point)

        # Handle index-based return (Shapely < 2.0 or specific configurations)
        # vs object-based return (Shapely 2.0+)
        if isinstance(result, (int, np.integer)):
            p_nearest_vertex = ref_vertices_objs[result]
        else:
            p_nearest_vertex = result

        dist_vertex = p_nearest_vertex.distance(point)

    # 3. Determine the closest candidate across both features
    # Start with the line as the candidate, override if vertex is closer
    p_final_nearest = p_nearest_line
    if dist_vertex < dist_line:
        p_final_nearest = p_nearest_vertex

    # 4. Apply the snapping strategy rules
    return _snapped_point_by_snapstrategy(
        point, p_final_nearest, p_nearest_vertex, snap_strategy, tolerance
    )


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


def extract_points_lines_from_geometry(geometry: ShapelyGeometry) -> GeometryCollection:
    """
    Recursively extracts all Points and LineStrings from a Shapely geometry object.

    Polygons are converted to their exterior and interior rings as LineStrings.
    Multi-geometries and GeometryCollections are traversed completely.

    Args:
        geometry: A Shapely geometry object (e.g., Point, Polygon, GeometryCollection).

    Returns:
        A GeometryCollection containing all extracted Point and LineString objects.
        Returns an empty GeometryCollection if the input geometry is empty.
    """
    extracted_geometries: List[base.BaseGeometry] = []

    if geometry.is_empty:
        # Return an empty GeometryCollection for an empty input
        return GeometryCollection()

    if isinstance(geometry, Polygon):
        # Convert Polygon to its exterior and interior rings (LineStrings)
        extracted_geometries.append(LineString(geometry.exterior.coords))
        for interior in geometry.interiors:
            extracted_geometries.append(LineString(interior.coords))

    elif isinstance(geometry, LineString):
        # LineString is returned as is
        extracted_geometries.append(geometry)

    elif isinstance(
        geometry, (MultiPoint, MultiPolygon, MultiLineString, GeometryCollection)
    ):
        # Recursively process components of collections/multi-geometries
        for geom in geometry.geoms:
            # We use .geoms to flatten the result from the recursive call
            extracted_geometries.extend(extract_points_lines_from_geometry(geom).geoms)

    elif isinstance(geometry, Point):
        # Point is returned as is
        extracted_geometries.append(geometry)

    # Return the collected geometries as a single GeometryCollection
    return GeometryCollection(extracted_geometries)


def to_crs(crs_input):
    """
    Converts a CRS input in any common format to a pyproj CRS

    Parameters:
        crs_input (int | str): A CRS code or name in any format, such as an integer (e.g., 4326),
                               a string like "EPSG:31370", "urn:ogc:def:crs:EPSG::3857", or "WGS 84".

    Returns:
        pyproj CRS object
    """
    try:
        crs = pyproj.CRS.from_user_input(crs_input)
        return crs
    except Exception as e:
        raise ValueError(f"Error interpreting CRS: {e}")


def from_crs(crs, format="uri"):
    """
    Converts a pyproj.CRS object to a specific representation.

    Parameters:
        crs (pyproj.CRS): A CRS object from pyproj.
        format (str): The desired output format. Options:
            - "uri": Returns the full PROJ string representation.
            - "id": Returns the EPSG integer code if available.
            - "epsg": Returns the authority string in the format "EPSG:xxxx".

    Returns:
        str | int | None:
            - str for "uri" or "epsg"
            - int for "id"
            - None if the requested representation is not available.

    Raises:
        ValueError: If conversion fails.
    """

    try:
        if format == "uri":
            auth = crs.to_authority()
            # return f"urn:ogc:def:crs:{auth[0]}::{auth[1]}"
            return f"http://www.opengis.net/def/crs/{auth[0]}/0/{auth[1]}"
        elif format == "id":
            return crs.to_epsg()
        elif format == "epsg":
            auth = crs.to_authority()
            return str(auth[0]) + ":" + str(auth[1])
    except Exception as e:
        raise ValueError(f"Error converting CRS: {e}")

def get_pseudo_coords (geom1,geom2):
    coords1 = get_coordinates(geom1)
    coords2 = get_coordinates(geom2)
    # (Tuples are 'hashable', NumPy arrays not)
    set1 = set(map(tuple, coords1))
    set2 = set(map(tuple, coords2))
    return set1 - set2

def get_non_pseudo_coords (geom1,geom2):
    coords1 = get_coordinates(geom1)
    coords2 = get_coordinates(geom2)
    # (Tuples are 'hashable', NumPy arrays not)
    set1 = set(map(tuple, coords1))
    set2 = set(map(tuple, coords2))
    pseudo = set1 - set2
    return set1 - pseudo


def multilinestring_multipoint_from_reference_intersection(
    reference_intersection
):
    if isinstance(reference_intersection, (LineString, MultiLineString)):
        return to_multi(reference_intersection), MultiPoint()
    if isinstance(reference_intersection, (Point, MultiPoint)):
        return MultiLineString(), to_multi(reference_intersection)
    if isinstance(reference_intersection, GeometryCollection):
        points = []
        lines = []
        for geom in reference_intersection.geoms:

            if isinstance(geom, (Point, MultiPoint)):
                points.append(geom)
            elif isinstance(geom, (LineString, MultiLineString)):
                lines.append(geom)
            elif isinstance(geom, (Polygon, MultiPolygon)):
                lines.append(geom.boundary)
            else:
                TypeError("Geometrytype not valid at this stage")
        return to_multi(safe_unary_union(lines)), to_multi(safe_unary_union(points))

    if isinstance(reference_intersection, (Polygon, MultiPolygon)):
        return to_multi(reference_intersection.boundary), MultiPoint()

    raise TypeError(
        "Reference could not be interpreted"
    )


def get_thematic_points(input_geometry, reference_intersection):
    """
    returns a MultiPoint geometry with points on the thematic geometry that are used for consistency while creating connection_lines
    :param input_geometry:
    :param reference_intersection:
    :return: MultiPoint
    """
    geom_to_process_line = input_geometry
    if isinstance(geom_to_process_line, LinearRing):
        geom_to_process_line = LineString(geom_to_process_line.coords)

    # geom_to_process_segmentized = segmentize(
    #     geom_to_process_line, self.config.partial_snap_max_segment_length
    # ) #We skip the segmentation
    geom_to_process_segmentized = geom_to_process_line
    # Split the line at all intersection points with the MultiLineString
    splitter = safe_unary_union(reference_intersection)
    try:
        geom_to_process_splitted = split(geom_to_process_segmentized, splitter)
    except (GeometryTypeError, ValueError):
        geom_to_process_splitted = geom_to_process_segmentized
    thematic_points = MultiPoint(
        list(get_coordinates(geom_to_process_splitted))
    )
    return thematic_points


def build_custom_network(
    input_geometry,
    theme_multiline,
    reference,
    reference_intersection,
    relevant_distance,
    snap_strategy=SnapStrategy.NO_PREFERENCE,
    gap_threshold=0.1,
    snap_dist=0.01,
):
    """
    Build a topological network from multiple linestrings and point geometries.

    This function integrates theme lines and reference lines into a graph,
    inserts pseudo-nodes (theme points) into the reference network, and
    executes a series of optimization steps to ensure connectivity and
    cycle formation.

    Parameters
    ----------
    input_geometry : shapely.geometry.MultiLineString
        The primary linestrings
    theme_multiline : shapely.geometry.MultiLineString
        The part of the input_geometry (linestrings) representing the thematic network layer.
    reference : shapely.geometry.MultiLineString
        The original reference linestrings used for snapping and connectivity.
    reference_intersection : shapely.geometry.MultiPoint
        The reference points used for snapping and connectivity.
    relevant_distance : float
        A base distance used to calculate interconnect and connectivity thresholds.
    snap_strategy: SnapStrategy
        Decides which connection lines between vertices are added.
    gap_threshold : float, optional
        The distance threshold for closing small gaps, by default 0.1.

    Returns
    -------
    nx.Graph
        A topologically connected graph with optimized edges and nodes.
    """
    G = nx.Graph()

    # TAG OVERVIEW
    # LINES:
    # ref_lines
    # theme_lines
    # interconnect
    # gap_closure
    # cycle_closure
    # component_interconnect

    # POINTS:
    # ref_points

    # VERTICES
    # theme_vertex
    # pseudo_theme_vertex
    # ref_vertex
    # pseudo_ref_vertex
    # pseudo_theme_ref_vertex

    ref_multiline, ref_points = (
        multilinestring_multipoint_from_reference_intersection(
            reference_intersection
        )
    )
    theme_points = get_thematic_points(input_geometry, reference_intersection)

    # 1. Load Theme & reference Lines
    pseudo_ref_coords = get_pseudo_coords(reference_intersection, reference)
    pseudo_theme_coords = get_pseudo_coords(theme_multiline, input_geometry)

    _multilinestring_to_edges(G, theme_multiline,pseudo_coords=pseudo_theme_coords,node_tag="theme_vertex",pseudonode_tag= "pseudo_theme_vertex", edge_tag = "theme_lines")
    srtree_ref_lines,edge_mapping_ref_lines =_multilinestring_to_edges(G, ref_multiline,pseudo_coords=pseudo_ref_coords,node_tag="ref_vertex",pseudonode_tag= "pseudo_ref_vertex", edge_tag = "ref_lines")

    # 2. Load Reference Points
    for pt in ref_points.geoms:
        p_coord = pt.coords[0]
        G.add_node(p_coord, tag="ref_points")

    # 3. Add Pseudo-nodes (theme_points) onto ref_lines
    if not ref_multiline is None and not ref_multiline.is_empty:
        for point in theme_points.geoms:
            if G.has_node(point.coords[0]):
                continue
            nearest_idx = srtree_ref_lines.nearest(
                point
            )
            # Find edge in NetworkX
            u, v = edge_mapping_ref_lines[nearest_idx]

            edge_data = G.get_edge_data(u, v)
            if not edge_data is None:
                line = edge_data["geometry"]
                p1, p2 = nearest_points(point, line)
            else:
                # This is possible when the original line is already split with a pseudonode
                edge_list = list(G.edges(u,data=True)) + list(G.edges(v,data=True))
                u,v,edge_data = find_closest_in_subset(point, edge_list)
                line = edge_data["geometry"]
                p1, p2 = nearest_points(point, line)
            if p2.distance(Point(u))>snap_dist and p2.distance(Point(v))>snap_dist:
                if p1.distance(p2) <= 1e-7:
                    _add_pseudonode(G, p2, u, v,tag_point="pseudo_theme_ref_vertex",tag_line="ref_lines")
                elif  p1.distance(p2) <= relevant_distance:
                    _add_pseudonode(G, p2, u, v,tag_point="pseudo_ref_vertex",tag_line="ref_lines")
    # 4. Optimization and Connectivity Pipeline
    # Close gaps and create initial theme-to-ref interconnections, keeping snapStrategy
    G = solve_all_network_gaps(
        G=G,
        snap_strategy=snap_strategy,
        gap_dist=gap_threshold,
        interconnect_dist=2*relevant_distance,
        snap_dist=snap_dist,
        merge_nodes=False,
    )

    # Ensure all isolated islands of reference lines are connected
    G = finalize_network_connectivity(
        G, interconnect_dist=2 * relevant_distance, edge_tags=["ref_lines"]
    )
    # remove pseudo_theme_vertices
    remove_pseodonodes (G,tag="pseudo_theme_vertex")

    # Force the creation of a cycle when input_geometry is a circle & the graph is currently a tree/forest
    if input_geometry.coords[0] == input_geometry.coords[-1]:
        G = ensure_cycle_exists(G, interconnect_dist=2 * relevant_distance)

    # # 5. Final Cleanup
    # composite_edges_removed = remove_composite_edges(G)
    # logging.debug(f"Removed composite edges: {len(composite_edges_removed)}")
    #export_to_geopackage(G,"g_out.gpkg")#TODO remove

    return G


def find_closest_in_subset(point, edge_list):
    """
    target_point_coords: (x, y) tuple
    edge_list: lijst van (u, v, data) triples uit G.edges(data=True)
    """
    best_dist = float("inf")
    closest_edge = None

    for u, v, data in edge_list:
        geom = data.get("geometry")
        if geom:
            # Bereken de kortste afstand van het punt tot deze specifieke lijn
            current_dist = point.distance(geom)

            if current_dist < best_dist:
                best_dist = current_dist
                closest_edge = (u, v, data)

    return closest_edge  # Geeft (u, v, data) terug van de winnaar


def remove_pseodonodes (G,tag):
    nodes_to_remove = [
        n for n, d in G.nodes(data=True) if d.get("tag") == tag
    ]

    for node in nodes_to_remove:
        # Get neighbors
        neighbors = list(G.neighbors(node))

        # Verbind alle buren met elkaar (bijv. bij een graad van 2 verbindt hij A met B)
        # Gebruik itertools als je nodes met een hogere graad dan 2 hebt


        for u, v in combinations(neighbors, 2):
            if not G.has_edge(u, v):
                G.add_edge(u, v, tag = "merged_lines",geometry=LineString([u, v]))

        # Remove old node
        G.remove_node(node)

def _add_pseudonode(G: Graph, p2, u, v,tag_point,tag_line):
    p_coord = p2.coords[0]
    # Remove the original edge and split it into two new edges at the point
    G.remove_edge(u, v)
    G.add_node(p_coord, tag=tag_point)
    G.add_edge(
        u, p_coord, tag=tag_line, geometry=LineString([u, p_coord])
    )
    G.add_edge(
        p_coord, v, tag=tag_line, geometry=LineString([p_coord, v])
    )


def _multilinestring_to_edges(G, multilinestring, node_tag,edge_tag,pseudo_coords,pseudonode_tag):
    multilinestring = to_multi(multilinestring)
    if not isinstance(multilinestring, MultiLineString):
        return None, None

    geoms_for_tree = []
    edge_mapping = []

    for line in multilinestring.geoms:
        coords = list(line.coords)
        for i in range(len(coords) - 1):
            u = coords[i]
            v = coords[i + 1]
            geom = LineString([u, v])

            # Voeg toe aan NetworkX
            G.add_edge(u, v, tag=edge_tag, geometry=geom)
            if not u in pseudo_coords:
                G.nodes[u]["tag"] = node_tag
            else:
                G.nodes[u]["tag"] = pseudonode_tag
            if not v in pseudo_coords:
                G.nodes[v]["tag"] = node_tag
            else:
                G.nodes[v]["tag"] = pseudonode_tag

            # Verzamel voor de STRtree
            geoms_for_tree.append(geom)
            edge_mapping.append((u, v))

    # Bouw de tree
    tree = STRtree(geoms_for_tree)

    # We retourneren de tree EN de mapping om de edges terug te kunnen vinden
    return tree, edge_mapping


def ensure_cycle_exists(G, interconnect_dist=5.0):
    """
    Ensure the graph contains at least one cycle by connecting the closest endpoints.

    If the graph is acyclic (a tree or forest), this function identifies all nodes
    with degree 1 and connects the two closest ones within the specified distance
    threshold to create a loop.

    Parameters
    ----------
    G : nx.Graph
        The network graph to evaluate.
    interconnect_dist : float, optional
        The maximum distance allowed to create a 'cycle_closure' edge between
        two endpoints, by default 5.0.

    Returns
    -------
    nx.Graph
        The graph, potentially with one additional 'cycle_closure' edge if a
        cycle was successfully formed.
    """
    # 1. Check if a cycle already exists
    try:
        cycles = nx.find_cycle(G)
        if cycles:
            logging.debug("Graph already contains a cycle. No action required.")
            return G
    except nx.NetworkXNoCycle:
        logging.debug("No cycle found. Searching for a closing link...")

    # 2. Identify all endpoints (nodes with degree 1)
    endpoints = [n for n, deg in G.degree() if deg == 1]

    if len(endpoints) < 2:
        logging.debug("Not enough endpoints to close a cycle.")
        return G

    # 3. Find the two closest endpoints using STRtree
    endpoint_pts = [Point(n) for n in endpoints]
    tree = STRtree(endpoint_pts)

    best_dist = float("inf")
    best_pair = None

    for i, pt_a in enumerate(endpoint_pts):
        # Search for candidates within the interconnect_dist buffer
        # Note: index 0 is often the node itself, handled by the 'i == idx' check
        indices = tree.query(pt_a.buffer(interconnect_dist))
        for idx in indices:
            if i == idx:
                continue

            pt_b = endpoint_pts[idx]
            dist = pt_a.distance(pt_b)

            if dist < best_dist:
                best_dist = dist
                best_pair = (endpoints[i], endpoints[idx])

    # 4. Add the "cycle-closing" edge
    if best_pair:
        u, v = best_pair
        G.add_edge(
            u, v, tag="cycle_closure", weight=best_dist, geometry=LineString([u, v])
        )
        logging.debug(f"Cycle closed between {u} and {v} (distance: {best_dist:.3f})")
    else:
        logging.debug("No endpoints found within the maximum distance to form a cycle.")

    return G


def solve_all_network_gaps(
    G, snap_strategy=SnapStrategy.NO_PREFERENCE, snap_dist=0.001, gap_dist=0.1, interconnect_dist=1.5, merge_nodes=False
):
    """
    Solve network gaps through universal snapping, gap closure, and targeted interconnects.

    Step 1 performs a universal gap closure across all layers of the graph. If enabled,
    it merges nodes within the snap distance. Step 2 creates functional connections
    specifically between 'theme_lines' endpoints and 'ref_lines'.

    Parameters
    ----------
    G : nx.Graph
        The input network graph containing nodes as coordinate tuples.
    snap_strategy : SnapStrategy
        The strategy to use to check which vertices should be connected
    snap_dist : float, optional
        Maximum distance for merging nodes (snapping), by default 0.001.
    gap_dist : float, optional
        Maximum distance for creating 'gap_closure' edges, by default 0.1.
    interconnect_dist : float, optional
        Maximum distance for 'interconnect' edges between theme and reference lines,
        by default 2.0.
    merge_nodes : bool, optional
        If True, nodes within `snap_dist` will be merged into one. If False,
        a 'gap_closure' edge is created instead, by default False.

    Returns
    -------
    nx.Graph
        The optimized graph with closed gaps and new interconnections.
    """
    # --- STEP 1: Universal Gap Closure ---
    all_nodes = list(G.nodes())
    if not all_nodes:
        return G

    all_points = [Point(n) for n in all_nodes]
    global_tree = STRtree(all_points)

    nodes_to_relabel = {}
    endpoints = [n for n, deg in G.degree() if deg <= 1]#endpoints or floating points

    for ep in endpoints:
        ep_pt = Point(ep)
        # Search for neighbors within gap_dist radius
        indices = global_tree.query(ep_pt.buffer(gap_dist))

        for idx in indices:
            candidate = all_nodes[idx]
            # Self-check: do not connect a node to itself
            if ep == candidate:
                continue

            dist = ep_pt.distance(all_points[idx])

            if dist <= gap_dist:
                if merge_nodes and dist <= snap_dist:
                    nodes_to_relabel[ep] = candidate
                    break
                else:
                    # Only add an edge if it does not form a self-loop
                    if ep != candidate and not G.has_edge(ep, candidate):
                        G.add_edge(
                            ep,
                            candidate,
                            tag="gap_closure",
                            weight=dist,
                            geometry=LineString([ep, candidate]),
                        )
                    break

    # Perform node merging (relabeling)
    if merge_nodes and nodes_to_relabel:
        # nx.relabel_nodes moves the connections.
        # If an edge existed between ep and candidate, it becomes a self-loop.
        nx.relabel_nodes(G, nodes_to_relabel, copy=False)

        # Explicitly remove all self-loops created by merging
        loops = list(nx.selfloop_edges(G))
        G.remove_edges_from(loops)

        if loops:
            logging.debug(f"Cleaned: {len(loops)} self-loops removed after merging.")

    #STEP 2 make interconnection edges
    theme_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get("tag") == "theme_lines"]
    ref_edges = [(u, v) for u, v, d in G.edges(data=True) if d.get("tag") == "ref_lines"]

    if ref_edges:
        theme_sub = G.edge_subgraph(theme_edges)
        ref_sub = G.edge_subgraph(ref_edges)

        theme_endpoints = [n for n, deg in theme_sub.degree() if deg == 1]

        nodes = list(ref_sub.nodes())
        trees = {}

        all_pts = [Point(n) for n in nodes]
        trees["all"] = {"pts": all_pts, "tree": STRtree(all_pts), "nodes": nodes}

        if snap_strategy != SnapStrategy.NO_PREFERENCE:
            v_nodes = [n for n in nodes if not G.nodes[n].get("tag", "").startswith("pseudo")]
            v_pts = [Point(n) for n in v_nodes]
            trees["vertex"] = {"pts": v_pts, "tree": STRtree(v_pts), "nodes": v_nodes}

        for tep in theme_endpoints:
            tep_pt = Point(tep)
            connection_made = False

            # Order of strategies to test
            strategies = ["vertex", "all"] if snap_strategy != SnapStrategy.NO_PREFERENCE else ["all"]

            for s_key in strategies:
                if s_key not in trees: continue

                s = trees[s_key]
                n_idx = s["tree"].nearest(tep_pt)
                if n_idx is None:
                    continue
                target_coord = s["nodes"][n_idx]
                dist = tep_pt.distance(s["pts"][n_idx])

                # Validation
                if tep != target_coord and gap_dist < dist <= interconnect_dist:
                    if not G.has_edge(tep, target_coord):
                        G.add_edge(
                            tep,
                            target_coord,
                            tag="interconnect",
                            weight=dist,
                            geometry=LineString([tep, target_coord]),
                        )
                        connection_made = True

                if connection_made:
                    break
    return G

def finalize_network_connectivity(
    G, interconnect_dist=2.0, edge_tags=["theme_lines", "ref_lines"]
):
    """
    Connect isolated network islands of reference lines if the graph is not fully connected.

    This function identifies disjoint components in the graph and attempts to bridge
    them by finding the nearest reference nodes between different components,
    provided they are within the specified distance threshold.

    Parameters
    ----------
    G : nx.Graph
        The network graph containing 'ref_lines' and 'theme_lines'.
    interconnect_dist : float, optional
        The maximum distance allowed to create a bridging connection between
        components, by default 2.0.

    Returns
    -------
    nx.Graph
        The graph with additional 'ref_interconnect' edges where gaps were bridged.
    """
    # 1. Retrieve all connected components
    components = list(nx.connected_components(G))
    if len(components) <= 1:
        return G  # Everything is already connected, no action needed

    # 2. Collect all potential reference nodes per component
    component_data = []
    for i, comp in enumerate(components):
        # Filter nodes that belong to a 'ref_lines' edge within this component
        ref_in_comp = [
            n
            for n in comp
            if any(d.get("tag") in edge_tags for _, _, d in G.edges(n, data=True))
        ]

        if ref_in_comp:
            comp_points = [Point(n) for n in ref_in_comp]
            component_data.append(
                {
                    "id": i,
                    "nodes": ref_in_comp,
                    "points": comp_points,
                    "tree": STRtree(comp_points),
                }
            )

    # 3. Attempt to connect the components
    for i in range(len(component_data)):
        comp_a = component_data[i]

        # Search for the nearest node in ALL other components
        best_dist = float("inf")
        best_connection = None

        for j in range(len(component_data)):
            if i == j:
                continue
            comp_b = component_data[j]

            # For each node in A, find the nearest neighbor in B using STRtree
            for idx_a, pt_a in enumerate(comp_a["points"]):
                idx_b = comp_b["tree"].nearest(pt_a)
                pt_b = comp_b["points"][idx_b]
                dist = pt_a.distance(pt_b)

                if dist < best_dist and dist <= interconnect_dist:
                    best_dist = dist
                    best_connection = (comp_a["nodes"][idx_a], comp_b["nodes"][idx_b])

        # 4. Add the interconnecting edge
        if best_connection:
            u, v = best_connection
            if not G.has_edge(u, v):
                G.add_edge(
                    u,
                    v,
                    tag="component_interconnect",
                    weight=best_dist,
                    geometry=LineString([u, v]),
                )
                logging.debug(
                    f"Component {i} connected to another component (distance: {best_dist:.3f})"
                )

    return G
