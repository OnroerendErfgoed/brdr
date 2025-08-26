from shapely import from_wkt
from shapely.ops import snap

from brdr.geometry_utils import safe_unary_union

# Dummy MultiLineString with small gaps
lines = from_wkt ("MULTILINESTRING ((171756.58074843569193035 171711.97158299898728728, 171769.32496534523670562 171701.38188981029088609, 171833.40245241450611502 171698.68555990885943174),(171791.03517037443816662 171700.68460767416399904, 171791.003030642750673 171700.4725628899759613, 171791.00269666060921736 171700.46974098280770704))")

print (lines.wkt)
# Step 1: Snap lines to themselves to close small gaps
snapped = snap(lines, lines, tolerance=0.001)
print (snapped.wkt)
unioned = safe_unary_union(snapped)
print (unioned.wkt)
