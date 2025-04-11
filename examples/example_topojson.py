
import json

import topojson
from shapely.geometry.linestring import LineString

from brdr.aligner import Aligner
from brdr.enums import GRBType
from brdr.geometry_utils import remove_shortest_and_merge
from brdr.grb import GRBActualLoader
from brdr.loader import DictLoader

with open("input/topo_parcels.geojson", 'r') as f:
    data = json.load(f)

assert data['type'] == 'FeatureCollection'
topo = topojson.Topology(data,prequantize=False)
print (topo.to_geojson())
print (topo.to_json())

new_arcs =[]
arc_id = 0
arc_dict ={}
for  arc in topo.output['arcs']:
    print(arc_id)
    print (arc)
    linestring = LineString(arc)
    print (linestring.wkt)
    print(linestring.geom_type)
    arc_dict[arc_id] = linestring
    arc_id = arc_id +1

aligner = Aligner(crs="EPSG:31370")
aligner.load_thematic_data(DictLoader(arc_dict))
aligner.load_reference_data(
    GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
)
relevant_distance = 0.4
process_result = aligner.process(
    relevant_distance=relevant_distance,
)
for arc_id in process_result.keys():
    print (str(arc_id))
    result_line = process_result[arc_id][relevant_distance]["result"]

    try:
        linestring = remove_shortest_and_merge(result_line)
        if linestring.geom_type=="MultiLineString":
            raise TypeError
        print("new_arc: " + linestring.wkt)
        new_arc = [list(coord) for coord in linestring.coords]
        new_arcs.append(new_arc)
    except:
        linestring = arc_dict[arc_id]
        print("old_arc: " + linestring.wkt)
        old_arc = [list(coord) for coord in linestring.coords]
        new_arcs.append(old_arc)


topo.output['arcs']=new_arcs

#todo uitzoeken
#Kunnen we multilinestrings toevoegen in arcs, of een boekhouding van arcs aanpassen per object?
#wat geeft het als je uitstekende linestrings in polygon samenvoegt: testen door topojson te manipuleren en dan om te zetten


# to visualize we use the (optional!) package Altair.
print(topo.to_geojson())