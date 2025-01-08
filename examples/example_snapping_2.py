from shapely import segmentize
from shapely.geometry import MultiPoint,Polygon, Point, LineString
from shapely.ops import nearest_points
from brdr.geometry_utils import get_coords_from_geometry
from shapely import segmentize
from shapely.geometry import MultiPoint, Polygon, Point, LineString
from shapely.ops import nearest_points


def trace_polygon(reference_poly, input_poly, tolerance):
    reference_poly_segmentized = segmentize(reference_poly, max_segment_length=1)
    # Ensure the vertices of the input polygon are on the boundary of the reference polygon
    reference_coords = MultiPoint(list(reference_poly_segmentized.exterior.coords))
    snapped_coords = []
    for coord in input_poly.exterior.coords:
        point = Point(coord)
        nearest_point = nearest_points(point, reference_coords)[1]
        if point.distance(nearest_point) <= tolerance:
            snapped_coords.append((nearest_point.x, nearest_point.y))
        else:
            snapped_coords.append((point.x, point.y))

    # Find the indices of the snapped points in the reference polygon boundary
    boundary_coords = list(reference_poly_segmentized.boundary.coords)
    snapped_indices = []
    for coord in snapped_coords:
        if (coord[0], coord[1]) in boundary_coords:
            snapped_indices.append(boundary_coords.index((coord[0], coord[1])))
        else:
            snapped_indices.append(None)

    # Create a new list of coordinates that includes the full parts of the reference polygon between the snapped vertices
    traced_coords = []
    for i in range(len(snapped_indices) - 1):
        start_index = snapped_indices[i]
        end_index = snapped_indices[i + 1]
        if start_index is not None and end_index is not None:
            if start_index < end_index:
                traced_coords.extend(boundary_coords[start_index: end_index + 1])
            else:
                traced_coords.extend(
                    boundary_coords[start_index:] + boundary_coords[: end_index + 1]
                )
        else:
            traced_coords.append(snapped_coords[i])
            traced_coords.append(snapped_coords[i + 1])

    # Create a new polygon from the traced coordinates
    traced_polygon = Polygon(traced_coords)

    return traced_polygon


def plot_polygons(reference_poly, input_poly, traced_poly):
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()

    # Plot reference polygon
    x1, y1 = reference_poly.exterior.xy
    ax.plot(x1, y1, 'b-', label='Reference Polygon')

    # Plot input polygon
    x2, y2 = input_poly.exterior.xy
    ax.plot(x2, y2, 'r--', label='Input Polygon')

    # Plot traced polygon
    x3, y3 = traced_poly.exterior.xy
    ax.plot(x3, y3, 'g--', label='Traced Polygon')

    ax.legend()
    plt.show()


# Example usage
reference_poly = Polygon([(0, 0), (4, 0), (4, 4), (0, 4), (0, 0)])
input_poly = Polygon([(1, 0.5), (3, 0.5), (3, 2), (1, 2), (1, 0.5)])
tolerance = 1

traced_poly = trace_polygon(reference_poly, input_poly, tolerance)

plot_polygons(reference_poly, input_poly, traced_poly)