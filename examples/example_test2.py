from brdr.grb import GRBActualLoader
# polygon
wkt="MULTIPOLYGON (((171795.71618631482124329 171817.88460136577486992, 171784.53532230854034424 171806.1688893586397171, 171746.73993028700351715 171841.20300138369202614, 171746.28380228579044342 171841.62578538432717323, 171767.35881029814481735 171856.89906539395451546, 171767.47471430152654648 171856.76376939192414284, 171798.38581831753253937 171820.68191336840391159, 171798.01820231974124908 171820.29669736698269844, 171795.71618631482124329 171817.88460136577486992)))"
from brdr.aligner import Aligner
from brdr.enums import OpenDomainStrategy, GRBType
from brdr.geometry_utils import geom_from_wkt
from brdr.loader import DictLoader, GeoJsonLoader

# CREATE AN ALIGNER
aligner = Aligner(
    crs="EPSG:31370",
)
# ADD A THEMATIC POLYGON TO THEMATIC DICTIONARY and LOAD into Aligner
thematic_dict = {"theme_id_1": geom_from_wkt(wkt)}
loader = DictLoader(thematic_dict)
aligner.load_thematic_data(loader)
# ADD A REFERENCE POLYGON TO REFERENCE DICTIONARY and LOAD into Aligner
reference_points ={
"type": "FeatureCollection",
"name": "points",
"crs": { "type": "name", "properties": { "name": "urn:ogc:def:crs:EPSG::31370" } },
"features": [
{ "type": "Feature", "properties": { "id": "2" }, "geometry": { "type": "Point", "coordinates": [ 171756.525063660374144, 171850.344575028371764 ] } },
{ "type": "Feature", "properties": { "id": "3" }, "geometry": { "type": "Point", "coordinates": [ 171766.127786894940073, 171857.041210968280211 ] } },
{ "type": "Feature", "properties": { "id": "4" }, "geometry": { "type": "Point", "coordinates": [ 171777.136171918798937, 171865.190890555531951 ] } },
{ "type": "Feature", "properties": { "id": "1" }, "geometry": { "type": "Point", "coordinates": [ 171745.722000021458371, 171841.942192198126577 ] } },
{ "type": "Feature", "properties": { "id": "6" }, "geometry": { "type": "Point", "coordinates": [ 171770.113516792538576, 171840.129037358070491 ] } }
]
}

reference_lines ={
"type": "FeatureCollection",
"name": "lines",
"crs": { "type": "name", "properties": { "name": "urn:ogc:def:crs:EPSG::31370" } },
"features": [
{ "type": "Feature", "properties": { "id": "1" }, "geometry": { "type": "LineString", "coordinates": [ [ 171741.11190000033821, 171839.010705479362514 ], [ 171751.689489041425986, 171847.237719177966937 ], [ 171762.267078082513763, 171855.699790410842979 ], [ 171772.727138356451178, 171862.986573972477345 ], [ 171783.539784931781469, 171870.8610013697471 ], [ 171794.94007534274715, 171877.79519862998859 ], [ 171801.874272603017744, 171884.259280821774155 ] ] } },
{ "type": "Feature", "properties": { "id": "2" }, "geometry": { "type": "LineString", "coordinates": [ [ 171761.437557435099734, 171790.265263877809048 ], [ 171748.806832151400158, 171775.059342177031795 ], [ 171759.598131422913866, 171766.107468917703955 ], [ 171736.053478466870729, 171738.393450334027875 ], [ 171722.809611179080093, 171721.715987823467003 ], [ 171702.208039842545986, 171713.867770171433222 ], [ 171657.080788343475433, 171750.165776812005788 ] ] } },
{ "type": "Feature", "properties": { "id": "3" }, "geometry": { "type": "LineString", "coordinates": [ [ 171657.238745108043076, 171750.592253888811683 ], [ 171683.091002234141342, 171729.768287099868758 ] ] } },
{ "type": "Feature", "properties": { "id": "4" }, "geometry": { "type": "LineString", "coordinates": [ [ 171683.456220801046584, 171729.549801494431449 ], [ 171656.795265419699717, 171750.054215320007643 ] ] } },
{ "type": "Feature", "properties": { "id": "6" }, "geometry": { "type": "LineString", "coordinates": [ [ 171758.284508655633545, 171787.105154872842832 ], [ 171747.304561027995078, 171775.410233911301475 ], [ 171758.641995322570438, 171765.655954856047174 ], [ 171749.347341981978389, 171756.25916246775887 ], [ 171724.833970534207765, 171776.074137721356237 ] ] } }
]
}


loader = GeoJsonLoader(_input=reference_points, id_property="id")
aligner.load_reference_data(loader)
# EXECUTE THE ALIGNMENT
relevant_distance = 3.2
process_result = aligner.process(
    relevant_distance=relevant_distance,
    od_strategy=OpenDomainStrategy.SNAP_ALL_SIDE,
    threshold_overlap_percentage=50,
)
# PRINT RESULTS IN WKT
print("result: " + process_result["theme_id_1"][relevant_distance]["result"].wkt)
print(
    "added area: "
    + process_result["theme_id_1"][relevant_distance]["result_diff_plus"].wkt
)
print(
    "removed area: "
    + process_result["theme_id_1"][relevant_distance]["result_diff_min"].wkt
)
