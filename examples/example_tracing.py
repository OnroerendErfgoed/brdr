import matplotlib.pyplot as plt

from brdr.enums import SnapStrategy
from brdr.geometry_utils import (
    snap_polygon_to_polygon,
    geom_from_wkt,
    polygon_to_multipolygon,
)

wkt_reference = "MultiPolygon (((170959.46072493179235607 174408.29865118899033405, 170973.48865984284202568 174389.69981613839627244, 170952.84080059177358635 174376.14473296594223939, 170939.2857174193195533 174374.25332601164700463, 170937.86716220359085128 174364.32343950160429813, 170922.73590656922897324 174354.39355299153248779, 170913.90934078249847516 174378.19375716644572094, 170959.46072493179235607 174408.29865118899033405)))"
geom_reference = geom_from_wkt(wkt_reference)
geom_reference = polygon_to_multipolygon(geom_reference)

wkt = "MultiPolygon (((170973.17342535045463592 174386.70508846075972542, 170928.56774467829382047 174355.96972545346943662, 170946.53611074411310256 174344.14843198910239153, 171004.06640560395317152 174369.36719137971522287, 170973.17342535045463592 174386.70508846075972542)))"
geom = geom_from_wkt(wkt)
geom = polygon_to_multipolygon(geom)

# Traceer rond de feature
traced_geom = snap_polygon_to_polygon(
    geometry=geom,
    reference=geom_reference,
    snap_strategy=SnapStrategy.PREFER_VERTICES,
    max_segment_length=5,
    tolerance=5,
)

# Visualisatie
fig, ax = plt.subplots()

# Plot de basispolygonen
for polygon in geom_reference.geoms:
    x, y = polygon.exterior.xy
    ax.plot(x, y, color="blue")

# Plot de traceerpolygonen
for polygon in geom.geoms:
    x, y = polygon.exterior.xy
    ax.plot(x, y, color="green")

# Plot het tracing resultaat
if traced_geom.geom_type == "Polygon":
    x, y = traced_geom.exterior.xy
    ax.plot(x, y, color="red")
elif traced_geom.geom_type == "MultiPolygon":
    for polygon in traced_geom.geoms:
        x, y = polygon.exterior.xy
        ax.plot(x, y, color="red")
if traced_geom.geom_type == "LineString":
    x, y = traced_geom.xy
    ax.plot(x, y, color="red")
elif traced_geom.geom_type == "MultiLineString":
    for line in traced_geom:
        x, y = line.xy
        ax.plot(x, y, color="red")
elif traced_geom.geom_type == "MultiPoint":
    for point in traced_geom.geoms:
        x = point.x
        y = point.y
        ax.plot(x, y, color="red")

ax.set_aspect("equal")
plt.show()
