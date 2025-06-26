from shapely import from_wkt, set_precision

wkt = "POINT (171761.43755743509973399 171790.2652638778090477)"
geom = from_wkt(wkt)
print(geom.wkt)
geom  = set_precision(geom, 0.01)
print (geom.wkt)
