#TESTSTRATEGY:
#featurecollection of geometries
#create topology
#adapt the arcs
#to_geojson


import topojson
import json

from shapely.geometry.linestring import LineString
from shapely.lib import line_merge, union
#from shapely.linear import line_merge
from shapely.ops import linemerge, unary_union
from shapely.set_operations import union_all

from brdr.aligner import Aligner
from brdr.enums import GRBType
from brdr.grb import GRBActualLoader
from brdr.loader import DictLoader

with open("input/one_simple.geojson", 'r') as f:
    data = json.load(f)

assert data['type'] == 'FeatureCollection'
topo = topojson.Topology(data,prequantize=False)
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
relevant_distance = 5
process_result = aligner.process(
    relevant_distance=relevant_distance,
)
for arc_id in process_result.keys():
    print (str(arc_id))
    result_line = process_result[arc_id][relevant_distance]["result"]

    linestring = linemerge(result_line)
    linestring_cleaned = []
    if linestring.geom_type=="MultiLineString":
        max_lines=0
        for line in linestring.geoms:
            print (line.wkt)
            new_arc = [list(coord) for coord in line.coords]
            new_arcs.append(new_arc)
            max_lines=max_lines +1
topo.output['arcs']=new_arcs
topo.output['objects']['data']['geometries'][0]['arcs']=[[list(range(0, max_lines))]]
print (topo)

#todo uitzoeken
#Kunnen we multilinestrings toevoegen in arcs, of een boekhouding van arcs aanpassen per object?
#wat geeft het als je uitstekende linestrings in polygon samenvoegt: testen door topojson te manipuleren en dan om te zetten


# to visualize we use the (optional!) package Altair.
print(topo.to_geojson())