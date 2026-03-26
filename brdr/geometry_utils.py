import json
import logging
import os
import tempfile
from io import BytesIO
from math import pi
from typing import Union, List, Tuple
from urllib.parse import urlencode

import geopandas as gpd
import numpy as np
import pyproj
import requests
from shapely import GEOSException, get_coordinates
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
from shapely import symmetric_difference
from shapely import to_wkt
from shapely import unary_union
from shapely import union
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
from shapely.ops import nearest_points
from shapely.ops import substring
from shapely.ops import unary_union
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
    if geometries is None:
        return GeometryCollection()
    if isinstance(geometries, BaseGeometry):
        return make_valid(geometries)
    if not isinstance(geometries, (list, tuple)):
        geometries = list(geometries)
    if len(geometries) == 0:
        return GeometryCollection()
    if len(geometries) == 1:
        return make_valid(geometries[0])
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
    angle_threshold_degrees=150.0,
):
    if geometry is None or geometry.is_empty or reference is None or reference.is_empty:
        return geometry
    geometry = to_multi(geometry, geomtype=None)
    if geometry.geom_type == "MultiPoint":
        result = _snap_point_to_reference(
            geometry,
            reference=reference,
            snap_strategy=snap_strategy,
            tolerance=tolerance,
            angle_threshold_degrees=angle_threshold_degrees,
        )
    elif geometry.geom_type == "MultiLineString":
        if max_segment_length > 0:
            geometry = _refine_multiline_near_reference_vertices(
                geometry,
                reference=reference,
                tolerance=tolerance,
                max_segment_length=max_segment_length,
            )
        result = _snap_line_to_reference(
            geometry,
            reference=reference,
            snap_strategy=snap_strategy,
            tolerance=tolerance,
            angle_threshold_degrees=angle_threshold_degrees,
        )
    elif geometry.geom_type == "MultiPolygon":
        if max_segment_length > 0:
            geometry = _refine_multipolygon_near_reference_vertices(
                geometry,
                reference=reference,
                tolerance=tolerance,
                max_segment_length=max_segment_length,
            )
        result = _snap_polygon_to_reference(
            geometry,
            reference=reference,
            snap_strategy=snap_strategy,
            tolerance=tolerance,
            angle_threshold_degrees=angle_threshold_degrees,
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
                angle_threshold_degrees=angle_threshold_degrees,
            )
            results.append(result)
        result = safe_unary_union(results)

    else:
        raise NotImplementedError(
            f"snapping for this type of geometry is not implemented: {str(geometry.geom_type)}"
        )
    return result


def _refine_multiline_near_reference_vertices(
    multiline: MultiLineString,
    reference: BaseGeometry,
    tolerance: float,
    max_segment_length: float,
) -> MultiLineString:
    if multiline is None or multiline.is_empty:
        return multiline
    ref_coords = list(get_coordinates(reference))
    if len(ref_coords) == 0:
        return multiline

    ref_vertex_points = [Point(c) for c in ref_coords]
    ref_vertex_tree = STRtree(ref_vertex_points) if ref_vertex_points else None

    refined_lines = []
    for line in multiline.geoms:
        refined_lines.append(
            _refine_line_near_reference_vertices(
                line,
                ref_vertex_points=ref_vertex_points,
                ref_vertex_tree=ref_vertex_tree,
                tolerance=tolerance,
                max_segment_length=max_segment_length,
            )
        )
    return MultiLineString(refined_lines)


def _refine_multipolygon_near_reference_vertices(
    multipolygon: MultiPolygon,
    reference: BaseGeometry,
    tolerance: float,
    max_segment_length: float,
) -> MultiPolygon:
    if multipolygon is None or multipolygon.is_empty:
        return multipolygon
    ref_coords = list(get_coordinates(reference))
    if len(ref_coords) == 0:
        return multipolygon

    ref_vertex_points = [Point(c) for c in ref_coords]
    ref_vertex_tree = STRtree(ref_vertex_points) if ref_vertex_points else None

    refined_polygons = []
    for polygon in multipolygon.geoms:
        # exterior
        exterior_line = LineString(polygon.exterior.coords)
        exterior_refined = _refine_line_near_reference_vertices(
            exterior_line,
            ref_vertex_points=ref_vertex_points,
            ref_vertex_tree=ref_vertex_tree,
            tolerance=tolerance,
            max_segment_length=max_segment_length,
        )
        exterior_coords = list(exterior_refined.coords)
        if exterior_coords[0] != exterior_coords[-1]:
            exterior_coords.append(exterior_coords[0])

        # interiors
        interior_rings = []
        for interior in polygon.interiors:
            interior_line = LineString(interior.coords)
            interior_refined = _refine_line_near_reference_vertices(
                interior_line,
                ref_vertex_points=ref_vertex_points,
                ref_vertex_tree=ref_vertex_tree,
                tolerance=tolerance,
                max_segment_length=max_segment_length,
            )
            interior_coords = list(interior_refined.coords)
            if interior_coords[0] != interior_coords[-1]:
                interior_coords.append(interior_coords[0])
            interior_rings.append(interior_coords)

        refined_polygon = make_valid(Polygon(exterior_coords, interior_rings))
        if refined_polygon is None or refined_polygon.is_empty:
            refined_polygon = polygon
        refined_polygons.append(refined_polygon)
    return to_multi(safe_unary_union(refined_polygons), geomtype="Polygon")


def _refine_line_near_reference_vertices(
    line: LineString,
    ref_vertex_points: List[Point],
    ref_vertex_tree: STRtree,
    tolerance: float,
    max_segment_length: float,
) -> LineString:
    if line is None or line.is_empty or line.length == 0:
        return line

    nearby_ref_points = []
    if ref_vertex_tree is not None:
        try:
            hits = ref_vertex_tree.query(line, predicate="dwithin", distance=tolerance)
            for h in hits:
                if isinstance(h, (int, np.integer)):
                    nearby_ref_points.append(ref_vertex_points[int(h)])
                else:
                    nearby_ref_points.append(h)
        except Exception:
            nearby_ref_points = [
                p for p in ref_vertex_points if line.distance(p) <= tolerance
            ]
    if not nearby_ref_points:
        return line

    # Keep original vertices and only add local refinement around projected ref vertices.
    sample_distances = [0.0, line.length]
    sample_distances.extend([line.project(Point(c)) for c in line.coords[1:-1]])

    half_step = max_segment_length / 2.0
    for p in nearby_ref_points:
        d = float(line.project(p))
        sample_distances.extend(
            [
                max(0.0, d - half_step),
                d,
                min(line.length, d + half_step),
            ]
        )

    sample_distances = sorted(set(round(float(d), 9) for d in sample_distances))
    new_coords = []
    for d in sample_distances:
        p = line.interpolate(d)
        coord = (float(p.x), float(p.y))
        if not new_coords or coord != new_coords[-1]:
            new_coords.append(coord)

    if len(new_coords) < 2:
        return line
    return LineString(new_coords)


def _snap_point_to_reference(
    geometry,
    reference,
    snap_strategy=SnapStrategy.NO_PREFERENCE,
    tolerance=1,
    angle_threshold_degrees=150.0,
):
    if geometry is None or geometry.is_empty or reference is None or reference.is_empty:
        return geometry

    if geometry.geom_type == "Point":
        geometry = MultiPoint([geometry])

    ref_border, ref_borders, ref_coords = _get_ref_objects(reference)
    ref_vertices_objs = [Point(c) for c in ref_coords]
    tree = STRtree(ref_vertices_objs) if ref_vertices_objs else None
    end_vertices, angle_vertices = _get_reference_vertex_priority_sets(
        ref_borders, angle_threshold_degrees=angle_threshold_degrees
    )

    if len(ref_coords) == 0:
        snap_strategy = SnapStrategy.NO_PREFERENCE

    points = []
    for geom in geometry.geoms:
        coords = list(geom.coords)
        coordinates = []
        for idx, coord in enumerate(coords):
            p = Point(coords[idx])
            p_snapped, bool_snapped, ref_vertices = _get_snapped_point(
                p,
                ref_border,
                tree,
                ref_vertices_objs,
                snap_strategy,
                tolerance,
                end_vertices=end_vertices,
                angle_vertices=angle_vertices,
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
    tolerance=1,
    angle_threshold_degrees=150.0,
):
    from brdr.graph_utils import longest_linestring_from_multilinestring

    # return snap(geometry,reference,tolerance)
    if geometry is None or geometry.is_empty or reference is None or reference.is_empty:
        return geometry

    if geometry.geom_type == "LineString":
        geometry = MultiLineString([geometry])

    ref_border, ref_borders, ref_coords = _get_ref_objects(reference)
    if len(ref_coords) == 0:
        snap_strategy = SnapStrategy.NO_PREFERENCE

    lines = []
    for geom in geometry.geoms:
        coords = list(geom.coords)
        coordinates = _get_snapped_coordinates(
            coords,
            ref_border,
            ref_borders,
            ref_coords,
            snap_strategy,
            tolerance,
            angle_threshold_degrees=angle_threshold_degrees,
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
    tolerance=1,
    correction_distance=0.01,
    angle_threshold_degrees=150.0,
):
    if geometry is None or geometry.is_empty or reference is None or reference.is_empty:
        return geometry
    geometry = to_multi(geometry, geomtype="Polygon")

    ref_border, ref_borders, ref_coords = _get_ref_objects(reference)
    if len(ref_coords) == 0:
        snap_strategy = SnapStrategy.NO_PREFERENCE
    polygons = []
    for geom in geometry.geoms:
        coords = list(geom.exterior.coords)
        coordinates = _get_snapped_coordinates(
            coords,
            ref_border,
            ref_borders,
            ref_coords,
            snap_strategy,
            tolerance,
            angle_threshold_degrees=angle_threshold_degrees,
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
    angle_threshold_degrees: float = 150.0,
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

    # 1. Lazy references for mixed-case handling
    ref_buffered_main = None
    ref_coords_buffered = None

    # Initialize Spatial Index (STRtree) for log(N) vertex lookups
    ref_vertices_objs = [Point(c) for c in ref_coords]
    tree = STRtree(ref_vertices_objs) if ref_vertices_objs else None
    end_vertices, angle_vertices = _get_reference_vertex_priority_sets(
        ref_borders, angle_threshold_degrees=angle_threshold_degrees
    )
    # Spatial index for nearest reference border lookup
    ref_lines_tree = STRtree(ref_borders) if ref_borders else None

    # 2. Cache the snapped start point for the first iteration
    p_start = Point(coords[0])
    p_start_snapped, bool_start_snapped, ref_vertices_start = _get_snapped_point(
        p_start,
        ref_border,
        tree,
        ref_vertices_objs,
        snap_strategy,
        tolerance,
        end_vertices=end_vertices,
        angle_vertices=angle_vertices,
    )

    for idx in range(1, len(coords)):
        p_end = Point(coords[idx])

        # Only snap the end point (the start point is carried over from the previous iteration)
        p_end_snapped, bool_end_snapped, ref_vertices_end = _get_snapped_point(
            p_end,
            ref_border,
            tree,
            ref_vertices_objs,
            snap_strategy,
            tolerance,
            end_vertices=end_vertices,
            angle_vertices=angle_vertices,
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
                    ref_lines_tree,
                    tolerance,
                )
            )
            coordinates.append(p_end_snapped.coords[0])

        else:
            # Mixed scenario: one point is snapped, the other is not
            line = LineString([p_start, p_end])

            # Select the appropriate pre-buffered geometry
            if ref_vertices_start or ref_vertices_end:
                if ref_coords_buffered is None:
                    ref_coords_multipoint = MultiPoint(ref_coords)
                    ref_coords_buffered = buffer_pos(ref_coords_multipoint, tolerance)
                curr_ref_buffered = ref_coords_buffered
            else:
                if ref_buffered_main is None:
                    ref_buffered_main = buffer_pos(ref_border, tolerance)
                curr_ref_buffered = ref_buffered_main
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
                        end_vertices=end_vertices,
                        angle_vertices=angle_vertices,
                    )

                    if first == p_start:
                        coordinates.extend(
                            _get_sublinestring_coordinates(
                                p_mid,
                                p_mid_snapped,
                                p_start,
                                p_start_snapped,
                                ref_borders,
                                ref_lines_tree,
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
                                ref_lines_tree,
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


def _as_coord_tuple(point: Point) -> Tuple[float, float]:
    return (float(point.x), float(point.y))


def _angle_between_vectors_degrees(vec_a: np.ndarray, vec_b: np.ndarray) -> float:
    """
    Compute the unsigned angle in degrees between two vectors.
    Returns 180.0 for degenerate vectors.
    """
    norm_a = float(np.linalg.norm(vec_a))
    norm_b = float(np.linalg.norm(vec_b))
    if norm_a == 0.0 or norm_b == 0.0:
        return 180.0
    cos_theta = float(np.dot(vec_a, vec_b) / (norm_a * norm_b))
    cos_theta = float(np.clip(cos_theta, -1.0, 1.0))
    return float(np.degrees(np.arccos(cos_theta)))


def _vertex_angle_degrees(a, b, c) -> float:
    ab = np.array([a[0] - b[0], a[1] - b[1]], dtype=float)
    cb = np.array([c[0] - b[0], c[1] - b[1]], dtype=float)
    return _angle_between_vectors_degrees(ab, cb)


def _get_reference_vertex_priority_sets(
    ref_borders: List[LineString], angle_threshold_degrees: float
) -> Tuple[set, set]:
    """
    Build coordinate sets for prioritized reference vertices:
    1) end vertices of open lines
    2) angle vertices with interior angle <= threshold
    """
    end_vertices = set()
    angle_vertices = set()

    for border in ref_borders:
        if border is None or border.is_empty or border.geom_type == "Point":
            continue
        coords = list(border.coords)
        if len(coords) < 2:
            continue

        is_closed = coords[0] == coords[-1]
        unique_coords = coords[:-1] if is_closed else coords
        n = len(unique_coords)
        if n < 2:
            continue

        if not is_closed:
            end_vertices.add((float(unique_coords[0][0]), float(unique_coords[0][1])))
            end_vertices.add((float(unique_coords[-1][0]), float(unique_coords[-1][1])))

        if n < 3:
            continue

        for i in range(n):
            if not is_closed and (i == 0 or i == n - 1):
                continue
            prev_i = (i - 1) % n
            next_i = (i + 1) % n
            angle = _vertex_angle_degrees(
                unique_coords[prev_i], unique_coords[i], unique_coords[next_i]
            )
            if angle <= angle_threshold_degrees:
                angle_vertices.add(
                    (float(unique_coords[i][0]), float(unique_coords[i][1]))
                )

    return end_vertices, angle_vertices


def _get_sublinestring_coordinates(
    p_end,
    p_end_snapped,
    p_start,
    p_start_snapped,
    ref_borders,
    ref_lines_tree,
    tolerance,
):
    coordinates = []
    if p_start_snapped == p_end_snapped:
        return coordinates
    reference_border, distance = _closest_line_fast(ref_borders, ref_lines_tree, p_end)
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


def _closest_line_fast(lines, lines_tree, point):
    if not lines:
        return None, -1
    if lines_tree is None:
        return closest_line(lines, point)
    nearest = lines_tree.nearest(point)
    if nearest is None:
        return None, -1
    if isinstance(nearest, (int, np.integer)):
        line = lines[int(nearest)]
    else:
        line = nearest
    return line, line.distance(point)


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


def _get_snapped_point(
    point: Point,
    ref_line: LineString,
    ref_coords_tree: STRtree,
    ref_vertices_objs: List[Point],
    snap_strategy: str,
    tolerance: float,
    end_vertices: set = None,
    angle_vertices: set = None,
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
        if snap_strategy == SnapStrategy.PREFER_VERTICES_ENDS_AND_ANGLES:
            candidate_vertices = []
            try:
                hits = ref_coords_tree.query(
                    point, predicate="dwithin", distance=tolerance
                )
                for h in hits:
                    if isinstance(h, (int, np.integer)):
                        candidate_vertices.append(ref_vertices_objs[int(h)])
                    else:
                        candidate_vertices.append(h)
            except Exception:
                pass

            if candidate_vertices:
                end_vertices = end_vertices or set()
                angle_vertices = angle_vertices or set()

                def _priority(candidate: Point):
                    coord = _as_coord_tuple(candidate)
                    if coord in end_vertices:
                        rank = 0
                    elif coord in angle_vertices:
                        rank = 1
                    else:
                        rank = 2
                    return rank, float(candidate.distance(point))

                p_nearest_vertex = min(candidate_vertices, key=_priority)
                dist_vertex = p_nearest_vertex.distance(point)
            else:
                result = ref_coords_tree.nearest(point)
                if isinstance(result, (int, np.integer)):
                    p_nearest_vertex = ref_vertices_objs[result]
                else:
                    p_nearest_vertex = result
                dist_vertex = p_nearest_vertex.distance(point)
        else:
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
    elif snap_strategy == SnapStrategy.PREFER_VERTICES_ENDS_AND_ANGLES:
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


def multilinestring_multipoint_from_reference_intersection(reference_intersection):
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

    raise TypeError("Reference could not be interpreted")







def gml_response_to_geojson(url, params, timeout=60):
    """
    Fetches a GML (Geography Markup Language) response from a URL,
    reads it using GeoPandas, and converts the result to
    GeoJSON format (as a Python dictionary).

    Args:
        url (str): The URL to fetch the GML data from (e.g., a WFS endpoint).
        params (dict): The parameters to be sent with the GET request
                       (e.g., service, version, request, typeName).

    Returns:
        dict: The geographical data in GeoJSON format.

    Raises:
        requests.exceptions.HTTPError: If the HTTP request fails.
        fiona.errors.DriverError: If the GML data cannot be read correctly.
    """
    # Fetch the GML response
    response = requests.get(url, params, timeout=timeout)
    response.raise_for_status()
    content = response.content

    read_attempts = [
        {"driver": "GML"},
        {"driver": "GML", "engine": "fiona"},
        {},
    ]

    last_exc = None

    # 1) Try reading directly from bytes
    for kwargs in read_attempts:
        try:
            gdf = gpd.read_file(BytesIO(content), **kwargs)
            return json.loads(gdf.to_json())
        except Exception as exc:
            last_exc = exc

    # 2) Fallback: write bytes to temporary .gml file and read explicitly
    tmp_path = None
    try:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".gml") as tmp:
            tmp.write(content)
            tmp_path = tmp.name

        for kwargs in read_attempts:
            try:
                gdf = gpd.read_file(tmp_path, **kwargs)
                return json.loads(gdf.to_json())
            except Exception as exc:
                last_exc = exc
    finally:
        if tmp_path and os.path.exists(tmp_path):
            os.remove(tmp_path)

    # 3) Fallback: let GDAL WFS driver fetch/parse the source directly
    try:
        query = urlencode(params or {}, doseq=True)
        full_url = f"{url}{'&' if '?' in url else '?'}{query}" if query else url

        wfs_read_attempts = [
            {"path": f"WFS:{full_url}", "kwargs": {"driver": "WFS"}},
            {"path": full_url, "kwargs": {"driver": "WFS"}},
            {"path": f"WFS:{full_url}", "kwargs": {"driver": "WFS", "engine": "fiona"}},
            {"path": full_url, "kwargs": {"driver": "WFS", "engine": "fiona"}},
            {"path": full_url, "kwargs": {}},
        ]

        for attempt in wfs_read_attempts:
            try:
                gdf = gpd.read_file(attempt["path"], **attempt["kwargs"])
                return json.loads(gdf.to_json())
            except Exception as exc:
                last_exc = exc
    except Exception as exc:
        last_exc = exc

    # All attempts failed
    raise last_exc


def total_vertex_distance(geom1, geom2, bidirectional=True) -> float:
    vertices1 = get_coordinates(geom1)
    total_distance_1 = 0.0
    len_vertices_1 = len(vertices1)
    for pt in vertices1:
        total_distance_1 += Point(pt).distance(geom2)
    total_distance = total_distance_1
    len_vertices = len_vertices_1

    if bidirectional:
        total_distance_2 = 0.0
        vertices2 = get_coordinates(geom2)
        len_vertices_2 = len(vertices2)
        for pt in vertices2:
            total_distance_2 += Point(pt).distance(geom1)
        total_distance = total_distance + total_distance_2
        len_vertices = len_vertices + len_vertices_2

    return total_distance / len_vertices
