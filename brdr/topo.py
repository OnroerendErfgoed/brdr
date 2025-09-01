import copy
import json

import topojson
from shapely import GeometryCollection, make_valid, LineString

from brdr.geometry_utils import (
    safe_unary_union,
    safe_difference,
    longest_linestring_from_multilinestring,
)
from brdr.utils import geojson_geometry_to_shapely


def dissolve_topo(
    dict_series,
    dict_thematic,
    dict_thematic_to_process,
    topo_thematic,
    relevant_distances,
):
    """
    Dissolves a processed dict_series of LineStrings (Arcs) into a dict_series of the original geometries
    :param dict_series:
    :param dict_thematic:
    :param dict_thematic_to_process:
    :param topo_thematic:
    :param relevant_distances:
    :return:
    """

    dict_series_topo = dict()
    for k, v in dict_thematic.items():
        dict_series_topo[k] = {}

    for relevant_distance in relevant_distances:
        for obj in topo_thematic.output["objects"]["data"]["geometries"]:
            key = obj["id"]
            topo = copy.deepcopy(topo_thematic)
            new_arcs = []
            for arc_id in dict_series.keys():
                try:
                    result_line = dict_series[arc_id][relevant_distance]["result"]

                    linestring = longest_linestring_from_multilinestring(result_line)
                    if linestring.geom_type == "MultiLineString":
                        raise TypeError
                    new_arc = [list(coord) for coord in linestring.coords]
                    new_arcs.append(new_arc)
                except:
                    linestring = dict_thematic_to_process[arc_id]
                    print("old_arc: " + linestring.wkt)
                    old_arc = [list(coord) for coord in linestring.coords]
                    new_arcs.append(old_arc)
            topo.output["arcs"] = new_arcs
            topo_geojson = topo.to_geojson()
            topo_geojson = json.loads(topo_geojson)
            result = GeometryCollection()
            for feature in topo_geojson["features"]:
                if feature["id"] == key:
                    result = geojson_geometry_to_shapely(feature["geometry"])
            result_diff_plus = make_valid(safe_difference(result, dict_thematic[key]))
            result_diff_min = make_valid(safe_difference(dict_thematic[key], result))
            result_diff = safe_unary_union([result_diff_plus, result_diff_min])
            dict_series_topo[key][relevant_distance] = {
                "result": result,
                "result_diff": result_diff,
                "result_diff_plus": result_diff_plus,
                "result_diff_min": result_diff_min,
                "result_relevant_intersection": GeometryCollection(),
                "result_relevant_diff": GeometryCollection(),
            }
    dict_series = dict_series_topo

    return dict_series


def generate_topo(dict_thematic_to_process):
    """
    Converts a dictionary (key-geometry) into a new dictionary (key-LineStrings) based on the arcs of a Topojson. It also returens the topojson
    :param dict_thematic_to_process: (key-geometry)
    :return: dict_thematic_to_process: (key-LineString (Arcs)) & Topojson
    """
    topo_thematic = topojson.Topology(dict_thematic_to_process, prequantize=False)
    print(topo_thematic.to_json())
    arc_id = 0
    arc_dict = {}
    for arc in topo_thematic.output["arcs"]:
        linestring = LineString(arc)
        arc_dict[arc_id] = linestring
        arc_id = arc_id + 1
    dict_thematic_to_process = arc_dict
    return dict_thematic_to_process, topo_thematic
