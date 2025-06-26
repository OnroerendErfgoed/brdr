from shapely import from_wkt, Polygon, LineString, make_valid
from shapely.geometry.polygon import LinearRing
from shapely.ops import linemerge

from brdr.geometry_utils import (
    buffer_pos,
    safe_intersection,
    safe_difference,
    safe_unary_union,
    to_multi,
    buffer_neg_pos,
)

thematic = from_wkt('MULTIPOLYGON (((69160.28090000152587891 225576.63089999929070473, 69153.27310000360012054 225567.70809999853372574, 69151.6039000004529953 225565.56210000067949295, 69149.9521000012755394 225563.40309999883174896, 69148.32169999927282333 225561.22749999910593033, 69146.70910000056028366 225559.03909999877214432, 69144.45549999922513962 225555.8977000005543232, 69142.24610000103712082 225552.72509999945759773, 69140.07190000265836716 225549.52809999883174896, 69137.91709999740123749 225546.31810000166296959, 69110.92010000348091125 225505.82009999826550484, 69086.86190000176429749 225469.40689999982714653, 69084.35729999840259552 225471.08109999820590019, 69081.85270000249147415 225472.75530000030994415, 69082.12110000103712082 225473.12209999933838844, 69115.85610000044107437 225523.86910000070929527, 69142.10209999978542328 225562.72109999880194664, 69150.0341000035405159 225573.30009999871253967, 69154.7450999990105629 225579.93409999832510948, 69155.04309999942779541 225580.38489999994635582, 69157.66210000216960907 225578.50789999961853027, 69160.28090000152587891 225576.63089999929070473)))')
reference = from_wkt('MULTILINESTRING ((69086.30299532413482666 225468.56099044904112816, 69083.89301132410764694 225470.3846704512834549, 69081.44961932301521301 225472.23369445279240608, 69073.84213131666183472 225460.82646244391798973, 69067.64501131325960159 225451.53398244082927704))')
relevant_distance = 2
buffer_distance = relevant_distance/2

thematic_exterior = thematic.geoms[0].exterior

thematic_buffered = buffer_pos(thematic_exterior,buffer_distance)
reference_buffered = buffer_pos(reference,buffer_distance)

overlap = safe_intersection(thematic_buffered,reference_buffered)
overlap_buffered = buffer_pos(overlap,buffer_distance)
reference_intersection=safe_intersection(reference,overlap_buffered)
thematic_difference=safe_difference(thematic_exterior,overlap_buffered)
line_list=[reference_intersection,thematic_difference]
#line_list=[thematic_difference,reference_intersection]
result = safe_unary_union(line_list)
result = linemerge(result)
result= to_multi(result)
coords = []
for r in result.geoms:
    coords.extend(r.coords)
# if coords[0] != coords[-1]:# linestring sluitend maken
#  coords.append(coords[0])
polygon = Polygon(coords)
line = LineString(coords)
linearring = LinearRing(coords)
polygon=buffer_neg_pos(polygon,0.01)
polygon=make_valid(polygon)
print (result.wkt)
print (polygon.wkt)
print (line.wkt)
print (linearring.wkt)

# for r
#
# coords = []
# for cols in csv_input:
#     coords.append((float(cols[0]), float(cols[1])))
#
# p = Polygon(coords)
#
# for feature in line_layer.getFeatures():
#  geom = feature.geometry()
#  if geom.isMultipart():
#  line_parts = geom.asMultiPolyline()
#  for part in line_parts:
#  points.extend(part)
#  else:
#  points.extend(geom.asPolyline())
#
#  # Ensure the polygon is closed by connecting the last point to the first if needed
#  if points[0] != points[-1]:
#  points.append(points[0])
#
#  # Create a polygon geometry
#  polygon_geom = QgsGeometry.fromPolygonXY([points])
#
