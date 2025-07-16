import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Voronoi
from shapely import from_wkt
from shapely.geometry import Polygon, LineString, MultiLineString
from shapely.ops import unary_union

# Step 1: Define a sample polygon
polygon = Polygon([(1, 1), (5, 1), (6, 3), (5, 5), (1, 5)])
polygon = from_wkt(
    "POLYGON ((171807.01653832942247391 171811.64178536087274551, 171806.50108232349157333 171811.20914536342024803, 171801.84252232313156128 171816.64696936681866646, 171798.38581831753253937 171820.68191336840391159, 171767.47471430152654648 171856.76376939192414284, 171767.35881029814481735 171856.89906539395451546, 171771.46146630495786667 171859.87231339514255524, 171777.14249030500650406 171863.9894334003329277, 171779.40213830769062042 171865.62706540152430534, 171783.0694023072719574 171860.27736939489841461, 171789.40066631138324738 171851.0414653904736042, 171790.78703431785106659 171849.01906538754701614, 171812.97858633100986481 171816.6465853676199913, 171809.78524232655763626 171813.96594536304473877, 171807.01653832942247391 171811.64178536087274551))"
)


# Step 2: Segmentize outline (interpolate points along the boundary)
def segmentize(polygon, segment_length=0.1):
    coords = list(polygon.exterior.coords)
    new_coords = []
    for i in range(len(coords) - 1):
        start = np.array(coords[i])
        end = np.array(coords[i + 1])
        dist = np.linalg.norm(end - start)
        num_segments = max(int(dist / segment_length), 1)
        for j in range(num_segments):
            point = start + (end - start) * j / num_segments
            new_coords.append(tuple(point))
    new_coords.append(coords[-1])
    return new_coords


segmented_points = segmentize(polygon)

# Step 3: Simplify if too many points
if len(segmented_points) > 500:
    polygon = polygon.simplify(0.5)
    segmented_points = segmentize(polygon)

# Step 4: Create Voronoi diagram
points = np.array(segmented_points)
vor = Voronoi(points)

# Step 5: Select Voronoi edges inside the polygon
lines = []
for vpair in vor.ridge_vertices:
    if -1 in vpair:
        continue
    p1 = vor.vertices[vpair[0]]
    p2 = vor.vertices[vpair[1]]
    line = LineString([p1, p2])
    if polygon.contains(line):
        lines.append(line)

# Step 6: Determine the best line (union of all valid Voronoi edges)
centerline = unary_union(lines)


# Step 7: Smooth the line (simple moving average for demonstration)
def smooth_line(line, window_size=3):
    if isinstance(line, LineString):
        coords = np.array(line.coords)
        if len(coords) < window_size:
            return line
        smoothed = []
        for i in range(len(coords)):
            start = max(0, i - window_size // 2)
            end = min(len(coords), i + window_size // 2 + 1)
            avg = np.mean(coords[start:end], axis=0)
            smoothed.append(tuple(avg))
        return LineString(smoothed)
    elif isinstance(line, MultiLineString):
        smoothed_lines = [smooth_line(ls, window_size) for ls in line.geoms]
        return MultiLineString(smoothed_lines)
    else:
        return line


smoothed_centerline = smooth_line(centerline)

# Plotting the result
x, y = polygon.exterior.xy
plt.plot(x, y, "black", label="Polygon")
for line in lines:
    x, y = line.xy
    plt.plot(x, y, "blue", alpha=0.5)
if isinstance(smoothed_centerline, LineString):
    x, y = smoothed_centerline.xy
    plt.plot(x, y, "red", linewidth=2, label="Smoothed Centerline")
elif isinstance(smoothed_centerline, MultiLineString):
    for line in smoothed_centerline.geoms:
        x, y = line.xy
        plt.plot(x, y, "red", linewidth=2)
plt.legend()
plt.axis("equal")
plt.title("Centerline Extraction from Polygon")
plt.show()
