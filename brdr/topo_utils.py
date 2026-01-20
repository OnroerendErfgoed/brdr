import copy
import json
from typing import Any, Dict

import topojson
from shapely import make_valid, LineString
from shapely.geometry import GeometryCollection
from shapely.geometry.base import BaseGeometry

from brdr.constants import REMARK_FIELD_NAME
from brdr.geometry_utils import (
    safe_unary_union,
    safe_difference,
    longest_linestring_from_multilinestring,
)
from brdr.typings import ProcessResult
from brdr.utils import geojson_geometry_to_shapely, union_process_result


def _dissolve_topo(
    thematic_id,
    process_results,
    input_geometry,
    thematic_geometries_to_process: Dict[Any, BaseGeometry],
    topo_thematic,
    relevant_distance,
) -> ProcessResult:
    """
    Dissolves a processed dictionary of LineStrings (Arcs) into a geometry collection
    of the original features based on TopoJSON.

    The function modifies the arcs in the provided TopoJSON object based on the
    processed geometries process results and then converts the resulting TopoJSON
    feature back into a Shapely geometry. It also calculates the symmetric
    difference between the resulting geometry and the input geometry.

    Args:
        thematic_id (str): The ID of the thematic object to dissolve in the TopoJSON.
        process_results (dict): A dictionary where keys are arc IDs and values contain
                            the processed result geometry (LineString or MultiLineString)
                            under the key relevant_distance -> "result".
        input_geometry (shapely.Geometry): The original input geometry for comparison.
        thematic_geometries_to_process (dict): A dictionary mapping arc IDs to their
                                         original Shapely LineString geometries.
        topo_thematic (TopoJSON object): The TopoJSON object containing the
                                         geometries and arcs to be processed.
        relevant_distance (float): The specific distance key used to retrieve the
                                   result from process_results.

    Returns:
        dict: A dictionary containing the resulting dissolved geometry and
              the various difference calculations:
              - "result": The dissolved Shapely geometry.
              - "result_diff": The symmetric difference (union of diff_plus and diff_min).
              - "result_diff_plus": Difference (result - input_geometry).
              - "result_diff_min": Difference (input_geometry - result).
              - "result_relevant_intersection": Empty GeometryCollection (placeholder).
              - "result_relevant_diff": Empty GeometryCollection (placeholder).
              - "properties": Dictionary for additional properties, including remarks.
    """

    # Create a deep copy of the TopoJSON object to avoid modifying the original
    topo = copy.deepcopy(topo_thematic)
    topo_geometries = topo.output["objects"]["data"]["geometries"]

    # Filter out geometries that do not match the target thematic_id
    # Note: Modifying a list while iterating requires careful handling (using a new list
    # or iterating over a copy/using list comprehension for removal).
    # The original code's removal logic is kept, but it is generally safer to iterate backwards
    # or build a new list. Assuming 'topo_geometries' is manageable:
    topo.output["objects"]["data"]["geometries"] = [
        obj for obj in topo_geometries if obj["id"] == thematic_id
    ]

    # The original geometry list was modified in place, so the loop structure is adjusted:
    if not topo.output["objects"]["data"]["geometries"]:
        # Handle case where thematic_id is not found, although the original code implies it should be present.
        pass

    # Process each arc and update the TopoJSON arcs list
    for obj in topo.output["objects"]["data"]["geometries"]:
        # Only process the target thematic object (only one should remain in the list)
        if obj["id"] != thematic_id:
            continue

        for arc_id in process_results.keys():
            try:
                # Get the processed line geometry
                result_line = process_results[arc_id][relevant_distance]["result"]

                # Extract the longest linestring if the result is a MultiLineString
                linestring = longest_linestring_from_multilinestring(result_line)

                if linestring.geom_type == "MultiLineString":
                    # Raise error if the resulting geometry is still a MultiLineString after processing
                    raise TypeError("Processed result is a MultiLineString.")

                # Convert the coordinates of the resulting linestring into the TopoJSON arc format
                new_arc = [list(coord) for coord in linestring.coords]

            except Exception:
                # If processing failed (e.g., TypeError or KeyError), use the original arc geometry
                linestring = thematic_geometries_to_process[arc_id]
                # Convert the coordinates of the original linestring into the TopoJSON arc format
                old_arc = [list(coord) for coord in linestring.coords]
                new_arc = old_arc

            # Update the arc in the TopoJSON object
            topo.output["arcs"][arc_id] = new_arc

    # Convert the modified TopoJSON object to GeoJSON string
    topo_geojson = topo.to_geojson()
    # Parse the GeoJSON string into a Python dictionary
    topo_geojson = json.loads(topo_geojson)
    result = GeometryCollection()

    # Extract the dissolved geometry corresponding to the thematic_id
    for feature in topo_geojson["features"]:
        if feature["id"] == thematic_id:
            # Convert the GeoJSON geometry structure to a Shapely geometry
            result = geojson_geometry_to_shapely(feature["geometry"])
            break  # Found the geometry, exit the loop

    # Calculate differences between the result and the input geometry
    # Difference (result - input_geometry)
    result_diff_plus = make_valid(safe_difference(result, input_geometry))
    # Difference (input_geometry - result)
    result_diff_min = make_valid(safe_difference(input_geometry, result))
    # Symmetric difference (union of the two differences)
    result_diff = safe_unary_union([result_diff_plus, result_diff_min])

    # Construct the final result dictionary
    process_result: ProcessResult = {
        "result": result,
        "result_diff": result_diff,
        "result_diff_plus": result_diff_plus,
        "result_diff_min": result_diff_min,
        "result_relevant_intersection": GeometryCollection(),  # Placeholder
        "result_relevant_diff": GeometryCollection(),  # Placeholder
        "properties": {REMARK_FIELD_NAME: []},
    }

    return union_process_result(process_result)


def _generate_topo(thematic_data):
    """
    Converts a dictionary of named Shapely geometries into a dictionary of
    LineStrings (Arcs) based on the generated TopoJSON structure.

    This process is essential for breaking down complex geometries into a
    shared set of simple edges (arcs) for further topological processing.

    Args:
        thematic_geometries_to_process (dict): A dictionary where keys are object IDs
                                         (e.g., thematic IDs) and values are
                                         Shapely geometries (e.g., Polygon, LineString).

    Returns:
        tuple: A tuple containing:
               - thematic_geometries_to_process (dict): A new dictionary where keys are
                 integer arc IDs and values are the corresponding Shapely LineStrings (Arcs).
               - topo_thematic (TopoJSON object): The generated TopoJSON object.
    """
    # Note: The original code commented out the deepcopy, suggesting it might
    # not be necessary if the input dictionary is intended to be overwritten
    # with arc geometries.

    # Generate the TopoJSON structure from the input geometries.
    # prequantize=False is used to prevent coordinate quantization.
    thematic_geometries: Dict[Any, BaseGeometry] = {}
    for key, feat in thematic_data.features.items():
        thematic_geometries[key] = feat.geometry
    topo_thematic = topojson.Topology(thematic_geometries, prequantize=False)

    # Print the resulting TopoJSON structure (for debugging/inspection)
    print(topo_thematic.to_json())

    arc_id = 0
    arc_dict = {}

    # Iterate through the generated arcs in the TopoJSON object
    for arc in topo_thematic.output["arcs"]:
        # Each arc is a list of coordinates; convert it into a Shapely LineString geometry
        linestring = LineString(arc)

        # Store the LineString in a dictionary, keyed by a sequential arc ID
        arc_dict[arc_id] = linestring
        arc_id = arc_id + 1

    # Replace the input dictionary with the new dictionary of arcs
    thematic_geometries_to_process = arc_dict

    # Return the dictionary of arcs and the TopoJSON object
    return thematic_geometries_to_process, topo_thematic


def _topojson_id_to_arcs(topojson):
    """
    Builds a dictionary mapping original feature IDs to their list of arc indices
    from a TopoJSON topology object.

    This function extracts the essential topological structure (how features are
    constructed from arcs) directly from the TopoJSON data.

    Args:
        topojson (TopoJSON object): The TopoJSON topology object, likely
                                    containing the structure under a '.output' attribute.

    Returns:
        dict: A dictionary where keys are the original feature IDs (as strings)
              and values are the list of arcs indices (the raw 'arcs' structure
              from the TopoJSON geometry).
    """
    dict_id_arcs = {}

    # Access the core data structure of the TopoJSON object
    topojson_data = topojson.output

    # TopoJSON features are contained under the "objects" key
    if "objects" not in topojson_data:
        print("Error: TopoJSON data does not contain an 'objects' section.")
        return dict_id_arcs

    # Iterate over all named geometry collections in 'objects' (e.g., 'data')
    for object_name, geo_collection in topojson_data["objects"].items():
        if "geometries" in geo_collection:
            # Iterate over each individual geometry feature within the collection
            for geometry in geo_collection["geometries"]:
                # Check if the geometry has both an 'id' and an 'arcs' definition
                if "id" in geometry and "arcs" in geometry:
                    original_id = geometry["id"]
                    # The 'arcs' data defines the geometry's shape using indices
                    # into the global TopoJSON arcs array. It can be an array of
                    # indices (LineString/Polygon) or an array of arrays of indices
                    # (MultiLineString/MultiPolygon). We store the full structure.
                    arcs_data = geometry["arcs"]

                    dict_id_arcs[original_id] = arcs_data

    return dict_id_arcs
