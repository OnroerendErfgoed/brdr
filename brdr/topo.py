import copy
import json

import topojson
from shapely import GeometryCollection, make_valid, LineString

from brdr.constants import REMARK_FIELD_NAME
from brdr.geometry_utils import (
    safe_unary_union,
    safe_difference,
    longest_linestring_from_multilinestring,
)
from brdr.utils import geojson_geometry_to_shapely


def dissolve_topo(
        thematic_id,
    dict_series,
    dict_thematic_topo_geoms,
    dict_thematic_to_process,
    topo_thematic,
    relevant_distance,
):
    """
    Dissolves a processed dict_series of LineStrings (Arcs) into a dict_series of the original geometries
    :param dict_series:
    :param dict_thematic_topo_geoms:
    :param dict_thematic_to_process:
    :param topo_thematic:
    :param relevant_distances:
    :return:
    """

    topo = copy.deepcopy(topo_thematic)
    topo_geometries = topo.output["objects"]["data"]["geometries"]
    for obj in topo_geometries:
        if obj["id"] != thematic_id:
            topo_geometries.remove(obj)
            continue
        for arc_id in dict_series.keys():
            try:
                result_line = dict_series[arc_id][relevant_distance]["result"]

                linestring = longest_linestring_from_multilinestring(result_line)
                if linestring.geom_type == "MultiLineString":
                    raise TypeError
                new_arc = [list(coord) for coord in linestring.coords]
            except:
                linestring = dict_thematic_to_process[arc_id]
                # print("old_arc: " + linestring.wkt)
                old_arc = [list(coord) for coord in linestring.coords]
                # new_arcs.append(old_arc)
                new_arc=old_arc
            topo.output["arcs"][arc_id]=new_arc
    topo_geojson = topo.to_geojson()
    topo_geojson = json.loads(topo_geojson)
    result = GeometryCollection()
    for feature in topo_geojson["features"]:
        if feature["id"] == thematic_id:
            result = geojson_geometry_to_shapely(feature["geometry"])
    result_diff_plus = make_valid(
        safe_difference(result, dict_thematic_topo_geoms[thematic_id])
    )
    result_diff_min = make_valid(
        safe_difference(dict_thematic_topo_geoms[thematic_id], result)
    )
    result_diff = safe_unary_union([result_diff_plus, result_diff_min])
    process_result = {
        "result": result,
        "result_diff": result_diff,
        "result_diff_plus": result_diff_plus,
        "result_diff_min": result_diff_min,
        "result_relevant_intersection": GeometryCollection(),
        "result_relevant_diff": GeometryCollection(),
        "properties": {REMARK_FIELD_NAME: []},
    }

    return process_result


def generate_topo(dict_thematic_to_process):
    """
    Converts a dictionary (key-geometry) into a new dictionary (key-LineStrings) based on the arcs of a Topojson. It also returns the topojson
    :param dict_thematic_to_process: (key-geometry)
    :return: dict_thematic_to_process: (key-LineString (Arcs)) & Topojson
    """
    dict_thematic_topo_geoms = copy.deepcopy(dict_thematic_to_process)
    topo_thematic = topojson.Topology(dict_thematic_to_process, prequantize=False)
    print(topo_thematic.to_json())
    arc_id = 0
    arc_dict = {}
    for arc in topo_thematic.output["arcs"]:
        linestring = LineString(arc)
        arc_dict[arc_id] = linestring
        arc_id = arc_id + 1
    dict_thematic_to_process = arc_dict
    return dict_thematic_to_process, topo_thematic, dict_thematic_topo_geoms


def topojson_id_to_arcs(topojson):
    """
    Bouwt een dictionary op van originele ID's naar de lijst van arcs-indices
    vanuit een TopoJSON-topologie.

    Args:
        topojson (dict): De TopoJSON-topologie geladen als een Python dictionary.

    Returns:
        dict: Een dictionary waarin de sleutels de originele ID's zijn (als string)
              en de waarden de lijst van arcs-indices.
    """
    dict_id_arcs = {}

    topojson_data=topojson.output
    # TopoJSON-data bevat de features onder de "objects" sleutel
    if "objects" not in topojson_data:
        print("Fout: De TopoJSON-data bevat geen 'objects' sectie.")
        return dict_id_arcs

    # Itereer over alle named geometry collections in 'objects'
    for object_name, geo_collection in topojson_data["objects"].items():
        if "geometries" in geo_collection:
            # Itereer over elke individuele geometrie
            for geometry in geo_collection["geometries"]:
                # We controleren of de geometrie een 'id' en 'arcs' heeft
                if "id" in geometry and "arcs" in geometry:
                    # De ID kan een string of een nummer zijn, we converteren naar string
                    original_id = str(geometry["id"])

                    # De arcs kunnen een array van arrays zijn (voor MultiPolygon/MultiLineString),
                    # of een array van indices (voor Polygon/LineString).
                    # We slaan de volledige 'arcs' structuur op zoals deze is.
                    arcs_data = geometry["arcs"]

                    dict_id_arcs[original_id] = arcs_data

    return dict_id_arcs
