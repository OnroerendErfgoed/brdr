import logging

import numpy as np
from shapely import GEOSException, from_wkt, to_wkt, STRtree
from shapely import GeometryCollection
from shapely import Polygon
from shapely import buffer
from shapely import difference
from shapely import intersection
from shapely import is_empty
from shapely import make_valid
from shapely import symmetric_difference
from shapely import unary_union
from shapely import union
from shapely.geometry import shape
from shapely.geometry.base import BaseGeometry
from shapely.prepared import prep

from brdr.constants import MITRE_LIMIT
from brdr.constants import QUAD_SEGMENTS
from brdr.constants import THRESHOLD_EXCLUSION_AREA
from brdr.constants import THRESHOLD_EXCLUSION_PERCENTAGE
from brdr.utils import get_collection


def buffer_neg_pos(geometry, buffer_value):
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
            quad_segs=QUAD_SEGMENTS,
            join_style="mitre",
            mitre_limit=MITRE_LIMIT,
        ),
        buffer_value,
        quad_segs=QUAD_SEGMENTS,
        join_style="mitre",
        mitre_limit=MITRE_LIMIT,
    )


def buffer_neg(geometry, buffer_value):
    """
    Computes the negative buffer of a given geometric object.

    Args:
        geometry (shapely.geometry.base.BaseGeometry): The input geometric object.
        buffer_value (float): The negative buffer distance.

    Returns:
        shapely.geometry.base.BaseGeometry: The result of applying the negative buffer.

    Notes:
        -   The function uses the Shapely library for geometric operations.
        -   Parameters like `quad_segs`, `join_style`, and `mitre_limit` are used for the
            buffer.

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
        quad_segs=QUAD_SEGMENTS,
        join_style="mitre",
        mitre_limit=MITRE_LIMIT,
    )


def buffer_pos(geometry, buffer_value):
    """
    Computes the positive buffer of a given geometric object.

    Args:
        geometry (shapely.geometry.base.BaseGeometry): The input geometric object.
        buffer_value (float): The positive buffer distance.

    Returns:
        shapely.geometry.base.BaseGeometry: The result of applying the positive buffer.

    Notes:
        -   The function uses the Shapely library for geometric operations.
        -   Parameters like `quad_segs`, `join_style`, and `mitre_limit` are used for the
            buffer.

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
        quad_segs=QUAD_SEGMENTS,
        join_style="mitre",
        mitre_limit=MITRE_LIMIT,
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
        except Exception:
            logging.error("error: empty geometry returned")
            geom = Polygon()

    return geom


def safe_intersection(geom_a: BaseGeometry, geom_b: BaseGeometry) -> BaseGeometry:
    """
    Calculates the intersection of two geometries with error handling.

    This function attempts to compute the intersection between two Shapely geometry objects (`geom_a` and `geom_b`).
    It incorporates error handling to address potential exceptions that might arise due to topological inconsistencies
    in the geometries, such as non-noded intersections between linestrings.

    Args:
        geom_a (BaseGeometry): The first Shapely geometry object.
        geom_b (BaseGeometry): The second Shapely geometry object.

    Returns:
        BaseGeometry: The intersection geometry as a Shapely object. It might be an empty Polygon if an error occurs during processing.

    Logs:
        - If a `GEOSException` occurs:
            - A warning message is logged with the WKT representations of both geometries.
            - The function attempts to buffer both geometries by a small value (0.0000001) and then perform the intersection.
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
        except Exception:
            logging.error("error: empty geometry returned")
            geom = Polygon()

    return geom


def safe_difference(geom_a, geom_b):
    """
    Calculates the difference between two geometries with error handling.

    This function computes the difference between two Shapely geometry objects (`geom_a` and `geom_b`).
    It incorporates error handling to address potential exceptions that might arise due to topological inconsistencies
    in the geometries, similar to non-noded intersections between linestrings.

    Args:
        geom_a (BaseGeometry): The first Shapely geometry object.
        geom_b (BaseGeometry): The second Shapely geometry object to be subtracted from the first.

    Returns:
        BaseGeometry: The difference geometry as a Shapely object. It might be an empty Polygon if an error occurs during processing.

    Logs:
        - If a `GEOSException` occurs:
            - A warning message is logged with the WKT representations of both geometries.
            - The function attempts to buffer both geometries by a small value (0.0000001) and then perform the difference operation.
        - If any other exception occurs:
            - An error message is logged indicating that an empty geometry is returned.
    """
    # function to solve exceptional error: shapely.errors.GEOSException:
    # TopologyException: found non-noded intersection between LINESTRING
    # see: https://gis.stackexchange.com/questions/50399
    try:
        geom = difference(geom_a, geom_b)
    except GEOSException:
        print("difference_error")
        try:
            logging.warning(
                "difference_error for geoms:" + geom_a.wkt + " and " + geom_b.wkt
            )
            geom = difference(buffer(geom_a, 0.0000001), buffer(geom_b, 0.0000001))
        except Exception:
            print("error: empty geometry returned")
            logging.error("error: empty geometry returned")
            geom = Polygon()

    return geom


def safe_symmetric_difference(geom_a, geom_b):
    """
    Calculates the symmetrical difference between two geometries with error handling.

    This function computes the symmetrical difference between two Shapely geometry objects (`geom_a` and `geom_b`).
    It incorporates error handling to address potential exceptions that might arise due to topological inconsistencies
    in the geometries, similar to non-noded intersections between linestrings.

    Args:
        geom_a (BaseGeometry): The first Shapely geometry object.
        geom_b (BaseGeometry): The second Shapely geometry object to be subtracted from the first.

    Returns:
        BaseGeometry: The symmetrical difference geometry as a Shapely object. It might be an empty Polygon if an error occurs during processing.

    Logs:
        - If a `GEOSException` occurs:
            - A warning message is logged with the WKT representations of both geometries.
            - The function attempts to buffer both geometries by a small value (0.0000001) and then perform the difference operation.
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
        except Exception:
            logging.error("error: empty geometry returned")
            geom = Polygon()

    return geom


def grid_bounds(geom: BaseGeometry, delta: float):
    """
    Divides a geometric area (specified by `geom`) into a grid of rectangular partitions.

    Args:
        geom (BaseGeometry): The geometric object representing the area to be partitioned.
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


def get_relevant_polygons_from_geom(geometry: BaseGeometry, buffer_distance: float):
    """
    Get only the relevant parts (polygon) from a geometry.
    Points, Lines and Polygons smaller than relevant distance are excluded from the result
    """
    if not geometry or geometry.is_empty:
        # If the input geometry is empty or None, do nothing.
        return geometry
    else:
        geometry = make_valid(unary_union(geometry))
        # Create a GeometryCollection from the input geometry.
        geometry_collection = GeometryCollection(geometry)
        array = []
        for g in geometry_collection.geoms:
            # Ensure each sub-geometry is valid.
            g = make_valid(g)
            if str(g.geom_type) in ["Polygon", "MultiPolygon"]:
                relevant_geom = buffer_neg(g, buffer_distance)
                if relevant_geom is not None and not relevant_geom.is_empty:
                    array.append(g)
    return make_valid(unary_union(array))


def calculate_geom_by_intersection_and_reference(
    geom_intersection: BaseGeometry,
    geom_reference: BaseGeometry,
    is_openbaar_domein,
    buffer_distance,
    threshold_overlap_percentage,
    threshold_exclusion_percentage=THRESHOLD_EXCLUSION_PERCENTAGE,
    threshold_exclusion_area=THRESHOLD_EXCLUSION_AREA,
):
    """
    Calculates the geometry based on intersection and reference geometries.

    Args:
        geom_intersection (BaseGeometry): The intersection geometry.
        geom_reference (BaseGeometry): The reference geometry.
        is_openbaar_domein (bool): A flag indicating whether it's a public domain
            (area not covered with reference polygon).
        threshold_exclusion_percentage (int): The threshold exclusion percentage.
        threshold_exclusion_area (int): The threshold exclusion area.
        buffer_distance (float): The buffer distance.
        threshold_overlap_percentage (int): The threshold overlap percentage.

    Returns:
        tuple: A tuple containing the resulting geometries:

        *   geom: BaseGeometry or None: The resulting geometry or None if conditions
            are not met.
        *   geom_relevant_intersection: BaseGeometry or None: The relevant
            intersection.
        *   geom_relevant_difference: BaseGeometry or None: The relevant difference.

    Notes:
        -   If the reference geometry area is 0, the overlap is set to 100%.
        -   If the overlap is less than relevant_OVERLAP_PERCENTAGE or the
            intersection area is less than relevant_OVERLAP_AREA, None is returned.
        -   Otherwise, the relevant intersection and difference geometries are
            calculated.
        -   If both relevant intersection and difference are non-empty, the final
            geometry is obtained by applying safe intersection and buffering.
        -   If only relevant intersection is non-empty, the result is the reference
            geometry.
        -   If only relevant difference is non-empty, the result is None.
    """

    if geom_reference.area == 0:
        overlap = 100

    else:
        overlap = geom_intersection.area * 100 / geom_reference.area

    if (
        overlap < threshold_exclusion_percentage
        or geom_intersection.area < threshold_exclusion_area
    ):
        return Polygon(), Polygon(), Polygon()

    geom_difference = safe_difference(geom_reference, geom_intersection)
    geom_relevant_intersection = buffer_neg(geom_intersection, buffer_distance)
    geom_relevant_difference = buffer_neg(geom_difference, buffer_distance)
    if (
        not geom_relevant_intersection.is_empty
        and not geom_relevant_difference.is_empty
    ):
        # relevant intersection and relevant difference
        geom_x = safe_intersection(
            geom_reference,
            safe_difference(
                geom_reference,
                safe_intersection(
                    geom_difference,
                    buffer_neg_pos(geom_difference, buffer_distance),
                ),
            ),
        )
        geom = safe_intersection(
            geom_x,
            buffer_pos(
                buffer_neg_pos(geom_x, buffer_distance),
                buffer_distance,
            ),
        )
        # when calculating for OD, we create a 'virtual parcel'. When calculating this virtual parcel, it is buffered to take outer boundaries into account.
        # This results in a side-effect that there are extra non-logical parts included in the result. The function below tries to exclude these non-logical parts.
        # see eo_id 206363 with relevant distance=0.2m and SNAP_ALL_SIDE
        if is_openbaar_domein:
            geom = get_relevant_polygons_from_geom(geom, buffer_distance)
    elif not geom_relevant_intersection.is_empty and geom_relevant_difference.is_empty:
        geom = geom_reference
    elif geom_relevant_intersection.is_empty and not geom_relevant_difference.is_empty:
        geom = geom_relevant_intersection  # (=empty geometry)
    else:
        if is_openbaar_domein:
            geom = geom_relevant_intersection  # (=empty geometry)
        # geom = snap_geom_to_reference (geom_intersection, geom_reference, relevant_distance)
        elif threshold_overlap_percentage < 0:
            # if we take a value of -1, the original border will be used
            geom = geom_intersection
        elif overlap > threshold_overlap_percentage:
            geom = geom_reference
        else:
            geom = geom_relevant_intersection  # (=empty geometry)
    return geom, geom_relevant_intersection, geom_relevant_difference


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


def create_donut(geometry, distance):
    inner_geometry = buffer_neg(geometry, distance)
    return safe_difference(geometry, inner_geometry)


def create_dictionary_from_url(
    bounds_array,
    buffer_value,
    output_dictionary,
    crs,
    geometry,
    key,
    limit,
    name_reference_id,
    prt,
    url,
):
    geom = buffer(
        geometry,
        buffer_value,
        quad_segs=QUAD_SEGMENTS,
        join_style="mitre",
        mitre_limit=MITRE_LIMIT,
    )
    bounds_array.append(geom)
    if prt < 1:
        output_dictionary = _create_dictionary(
            output_dictionary, crs, geom, key, limit, name_reference_id, url
        )
    if prt > 0:
        geom = unary_union(bounds_array)
        grid = partition(geom, prt)
        for g in grid:
            output_dictionary = _create_dictionary(
                output_dictionary, crs, g, key, limit, name_reference_id, url
            )
    return output_dictionary


def _create_dictionary(input_dict, crs, geom, key, limit, name_reference_id, url):
    output_dict = {}
    output_dict.update(input_dict)
    bbox = str(geom.bounds).strip("()")
    url_with_bbox = (
        url
        + "&f=application%2Fgeo%2Bjson"
        + "&limit="
        + str(limit)
        + "&crs="
        + crs
        + "&bbox-crs="
        + crs
        + "&bbox="
        + bbox
    )
    logging.debug(key + "-->" + str(url_with_bbox))
    coll = _get_dict_from_url(url_with_bbox, name_reference_id, limit)
    output_dict.update(coll)
    return output_dict


def _get_dict_from_url(input_url, name_reference_id, limit):
    collection = get_collection(input_url, limit)
    dictionary = {}
    if "features" not in collection or len(collection["features"]) == 0:
        return dictionary
    for f in collection["features"]:
        key = str(f["properties"][name_reference_id])
        geom = shape(f["geometry"])
        if key not in collection:
            dictionary[key] = make_valid(geom)
        logging.debug(key + "-->" + str(geom))
    return dictionary


def features_by_geometric_operation(
    list_input_geometries, list_input_ids, list_geometries, predicate="intersects"
):
    thematic_tree = STRtree(list_input_geometries)
    thematic_items = np.array(list_input_ids)
    arr_indices = thematic_tree.query(list_geometries, predicate=predicate)
    thematic_intersections = list(set(thematic_items.take(arr_indices[1])))
    return thematic_intersections


def partition(geom, delta):
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
    prepared_geom = prep(geom)
    partitions = grid_bounds(geom, delta)
    filtered_grid = list(filter(prepared_geom.intersects, partitions))
    return filtered_grid
