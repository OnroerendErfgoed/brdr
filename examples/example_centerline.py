import matplotlib.pyplot as plt
from centerline.geometry import Centerline
from shapely import from_wkt
from shapely.geometry import LineString

# polygon = Polygon([[0, 0], [0, 4], [4, 4], [4, 0]])
polygon = from_wkt(
    "POLYGON ((171807.01653832942247391 171811.64178536087274551, 171806.50108232349157333 171811.20914536342024803, 171801.84252232313156128 171816.64696936681866646, 171798.38581831753253937 171820.68191336840391159, 171767.47471430152654648 171856.76376939192414284, 171767.35881029814481735 171856.89906539395451546, 171771.46146630495786667 171859.87231339514255524, 171777.14249030500650406 171863.9894334003329277, 171779.40213830769062042 171865.62706540152430534, 171783.0694023072719574 171860.27736939489841461, 171789.40066631138324738 171851.0414653904736042, 171790.78703431785106659 171849.01906538754701614, 171812.97858633100986481 171816.6465853676199913, 171809.78524232655763626 171813.96594536304473877, 171807.01653832942247391 171811.64178536087274551))"
)
attributes = {"id": 1, "name": "polygon", "valid": True}
centerline = Centerline(polygon, **attributes)
print(centerline.id == 1)
print(centerline.name)
print(centerline.geometry.geoms)
smoothed_centerline = centerline.geometry
lines = centerline.geometry.geoms

# Plotting the result
x, y = polygon.exterior.xy
plt.plot(x, y, "black", label="Polygon")
for line in lines:
    x, y = line.xy
    plt.plot(x, y, "blue", alpha=0.5)
if isinstance(smoothed_centerline, LineString):
    x, y = smoothed_centerline.xy
    plt.plot(x, y, "red", linewidth=2, label="Smoothed Centerline")
# elif isinstance(smoothed_centerline, MultiLineString):
#     for line in smoothed_centerline.geoms:
#         x, y = line.xy
#         plt.plot(x, y, 'red', linewidth=2)
plt.legend()
plt.axis("equal")
plt.title("Centerline Extraction from Polygon")
plt.show()
