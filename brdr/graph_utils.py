import logging
from itertools import combinations
from math import inf

import networkx as nx
import numpy as np
from networkx import Graph
from shapely import STRtree, get_coordinates, polygonize, unary_union
from shapely.errors import GeometryTypeError
from shapely.geometry import (
    GeometryCollection,
    LinearRing,
    LineString,
    MultiLineString,
    MultiPoint,
    MultiPolygon,
    Point,
    Polygon,
)
from shapely.ops import nearest_points, split

from brdr.enums import SnapStrategy
from brdr.geometry_utils import (
    _angle_between_vectors_degrees,
    safe_unary_union,
    total_vertex_distance,
    to_multi,
)


def _is_angle_graph_node(graph, node_id, angle_threshold_degrees=150.0):
    neighbors = list(graph.neighbors(node_id))
    if len(neighbors) < 2:
        return False
    center = np.array([float(node_id[0]), float(node_id[1])], dtype=float)
    vectors = []
    for nb in neighbors:
        vec = np.array([float(nb[0]), float(nb[1])], dtype=float) - center
        norm = np.linalg.norm(vec)
        if norm > 0:
            vectors.append(vec / norm)
    if len(vectors) < 2:
        return False
    min_angle = 180.0
    for i in range(len(vectors)):
        for j in range(i + 1, len(vectors)):
            angle = _angle_between_vectors_degrees(vectors[i], vectors[j])
            if angle < min_angle:
                min_angle = angle
    return min_angle <= angle_threshold_degrees


def nearest_node(point, nodes):
    return min(nodes, key=lambda n: Point(n).distance(point))


def get_pseudo_coords(geom1, geom2):
    coords1 = get_coordinates(geom1)
    coords2 = get_coordinates(geom2)
    set1 = set(map(tuple, coords1))
    set2 = set(map(tuple, coords2))
    return set1 - set2


def get_end_coords(geom1, geom2=None):
    geom1 = to_multi(geom1, geomtype=None)
    end_coords = set()

    if isinstance(geom1, LineString):
        coords = list(geom1.coords)
        if coords:
            end_coords.add(tuple(coords[0]))
            end_coords.add(tuple(coords[-1]))
    elif isinstance(geom1, MultiLineString):
        for line in geom1.geoms:
            coords = list(line.coords)
            if coords:
                end_coords.add(tuple(coords[0]))
                end_coords.add(tuple(coords[-1]))

    if geom2 is None:
        return end_coords

    coords2 = get_coordinates(geom2)
    set2 = set(map(tuple, coords2))
    return end_coords.intersection(set2)


def multilinestring_multipoint_from_reference_intersection(reference_intersection):
    if isinstance(reference_intersection, (LineString, MultiLineString)):
        return to_multi(reference_intersection), MultiPoint()
    if isinstance(reference_intersection, (Point, MultiPoint)):
        return MultiLineString(), to_multi(reference_intersection)
    if isinstance(reference_intersection, GeometryCollection):
        points = []
        lines = []
        for geom in reference_intersection.geoms:
            if isinstance(geom, (Point, MultiPoint)):
                points.append(geom)
            elif isinstance(geom, (LineString, MultiLineString)):
                lines.append(geom)
            elif isinstance(geom, (Polygon, MultiPolygon)):
                lines.append(geom.boundary)
            else:
                TypeError("Geometrytype not valid at this stage")
        return to_multi(safe_unary_union(lines)), to_multi(safe_unary_union(points))

    if isinstance(reference_intersection, (Polygon, MultiPolygon)):
        return to_multi(reference_intersection.boundary), MultiPoint()

    raise TypeError("Reference could not be interpreted")


def get_thematic_points(input_geometry, reference_intersection):
    geom_to_process_line = input_geometry
    if isinstance(geom_to_process_line, LinearRing):
        geom_to_process_line = LineString(geom_to_process_line.coords)
    geom_to_process_segmentized = geom_to_process_line
    splitter = safe_unary_union(reference_intersection)
    try:
        geom_to_process_splitted = split(geom_to_process_segmentized, splitter)
    except (GeometryTypeError, ValueError):
        geom_to_process_splitted = geom_to_process_segmentized
    thematic_points = MultiPoint(list(get_coordinates(geom_to_process_splitted)))
    return thematic_points


def _add_pseudonode(G: Graph, p2, u, v, tag_point, tag_line):
    p_coord = p2.coords[0]
    G.remove_edge(u, v)
    G.add_node(p_coord, tag=tag_point)
    geom = LineString([u, p_coord])
    G.add_edge(u, p_coord, tag=tag_line, geometry=geom, length=geom.length)
    geom = LineString([p_coord, v])
    G.add_edge(p_coord, v, tag=tag_line, geometry=geom, length=geom.length)


def find_closest_in_subset(target_point, edge_list):
    geoms = [data.get("geometry") for _, _, data in edge_list if data.get("geometry")]
    if not geoms:
        return None
    tree = STRtree(geoms)
    nearest_idx = tree.nearest(target_point)
    return edge_list[nearest_idx]


def add_pseudonodes_to_ref_line(
    G,
    edge_mapping_ref_lines,
    ref_multiline,
    relevant_distance: float,
    snap_dist: float,
    srtree_ref_lines,
    theme_points: MultiPoint,
):
    if not ref_multiline is None and not ref_multiline.is_empty:
        for point in theme_points.geoms:
            if G.has_node(point.coords[0]):
                continue
            nearest_idx = srtree_ref_lines.nearest(point)
            u, v = edge_mapping_ref_lines[nearest_idx]

            edge_data = G.get_edge_data(u, v)
            if not edge_data is None:
                line = edge_data["geometry"]
                p1, p2 = nearest_points(point, line)
            else:
                edge_list = list(G.edges(u, data=True)) + list(G.edges(v, data=True))
                u, v, edge_data = find_closest_in_subset(point, edge_list)
                line = edge_data["geometry"]
                p1, p2 = nearest_points(point, line)
            if p2.distance(Point(u)) > snap_dist and p2.distance(Point(v)) > snap_dist:
                if p1.distance(p2) <= 1e-7:
                    _add_pseudonode(
                        G,
                        p2,
                        u,
                        v,
                        tag_point="pseudo_ref_vertex_1",
                        tag_line="ref_lines",
                    )
                elif p1.distance(p2) <= relevant_distance:
                    _add_pseudonode(
                        G,
                        p2,
                        u,
                        v,
                        tag_point="pseudo_ref_vertex_2",
                        tag_line="ref_lines",
                    )
            elif (
                p2.distance(Point(u)) <= snap_dist
                and G.nodes[u]["tag"] == "pseudo_ref_vertex"
            ):
                G.nodes[u]["tag"] = "pseudo_ref_vertex_3"
            elif (
                p2.distance(Point(v)) <= snap_dist
                and G.nodes[v]["tag"] == "pseudo_ref_vertex"
            ):
                G.nodes[v]["tag"] = "pseudo_ref_vertex_3"
    return G


def remove_pseudonodes(G, tag, target_degrees=[2]):
    nodes_to_remove = [
        n
        for n, d in G.nodes(data=True)
        if d.get("tag") == tag and G.degree(n) in target_degrees
    ]

    for node in nodes_to_remove:
        neighbors = list(G.neighbors(node))
        for u, v in combinations(neighbors, 2):
            if not G.has_edge(u, v):
                pos_u = Point(u)
                pos_v = Point(v)
                new_line = LineString([pos_u, pos_v])
                G.add_edge(
                    u,
                    v,
                    tag="connection_" + tag + "_removed",
                    geometry=new_line,
                    length=new_line.length,
                )
        G.remove_node(node)
    return G


def _multilinestring_to_edges(
    G,
    multilinestring,
    node_tag,
    edge_tag,
    pseudo_coords,
    pseudonode_tag,
    end_coords=None,
    endnode_tag=None,
):
    multilinestring = to_multi(multilinestring)
    if not isinstance(multilinestring, MultiLineString):
        return None, None

    geoms_for_tree = []
    edge_mapping = []
    end_coords = set(end_coords or [])

    for line in multilinestring.geoms:
        coords = list(line.coords)
        for i in range(len(coords) - 1):
            u = coords[i]
            v = coords[i + 1]
            geom = LineString([u, v])

            G.add_edge(u, v, tag=edge_tag, geometry=geom, length=geom.length)
            if not u in pseudo_coords:
                G.nodes[u]["tag"] = node_tag
            else:
                G.nodes[u]["tag"] = pseudonode_tag
            if not v in pseudo_coords:
                G.nodes[v]["tag"] = node_tag
            else:
                G.nodes[v]["tag"] = pseudonode_tag
            if endnode_tag:
                if u in end_coords:
                    G.nodes[u][endnode_tag] = True
                if v in end_coords:
                    G.nodes[v][endnode_tag] = True

            geoms_for_tree.append(geom)
            edge_mapping.append((u, v))
    tree = STRtree(geoms_for_tree)
    return tree, edge_mapping


def connect_network_gaps(G, snap_dist=0.001, gap_dist=0.1, merge_nodes=False):
    all_nodes = list(G.nodes())
    if not all_nodes:
        return G

    all_points = [Point(n) for n in all_nodes]
    global_tree = STRtree(all_points)
    nodes_to_relabel = {}
    endpoints = [n for n, deg in G.degree() if deg <= 1]

    for ep in endpoints:
        ep_pt = Point(ep)
        indices = global_tree.query(ep_pt, predicate="dwithin", distance=gap_dist)

        for idx in indices:
            candidate = all_nodes[idx]
            if ep == candidate:
                continue

            dist = ep_pt.distance(all_points[idx])

            if merge_nodes and dist <= snap_dist:
                nodes_to_relabel[ep] = candidate
                break
            else:
                if ep != candidate and not G.has_edge(ep, candidate):
                    geom = LineString([ep, candidate])
                    G.add_edge(
                        ep,
                        candidate,
                        tag="gap_closure",
                        geometry=geom,
                        length=geom.length,
                    )
                break

    if merge_nodes and nodes_to_relabel:
        nx.relabel_nodes(G, nodes_to_relabel, copy=False)
        loops = list(nx.selfloop_edges(G))
        G.remove_edges_from(loops)
        if loops:
            logging.debug(f"Cleaned: {len(loops)} self-loops removed after merging.")
    return G


def connect_theme_and_reference(G, gap_dist=0.1, interconnect_dist=1.5):
    theme_edges = [
        (u, v) for u, v, d in G.edges(data=True) if d.get("tag") == "theme_lines"
    ]
    if not theme_edges:
        return G
    theme_sub = G.edge_subgraph(theme_edges)
    theme_endpoints = [n for n, deg in theme_sub.degree() if deg == 1]
    if not theme_endpoints:
        return G
    ref_node_ids = {
        node
        for u, v, d in G.edges(data=True)
        if d.get("tag") in ("ref_lines", "ref_points_connection")
        for node in (u, v)
    }
    extra_ref_nodes = [u for u, d in G.nodes(data=True) if d.get("tag") == "ref_points"]
    ref_node_ids.update(extra_ref_nodes)
    ref_sub_graph = nx.Graph()
    if ref_node_ids:
        ref_sub_graph.add_nodes_from((n, G.nodes[n]) for n in ref_node_ids if n in G)
    else:
        return G
    nodes = list(ref_sub_graph.nodes())
    all_pts = [Point(n) for n in nodes]
    all_tree = STRtree(all_pts)

    for tep in theme_endpoints:
        tep_pt = Point(tep)
        n_idx = all_tree.nearest(tep_pt)
        if n_idx is None:
            continue
        target_coord = nodes[n_idx]
        dist = tep_pt.distance(all_pts[n_idx])
        if tep != target_coord and gap_dist < dist <= interconnect_dist:
            if not G.has_edge(tep, target_coord):
                geom = LineString([tep, target_coord])
                G.add_edge(
                    tep,
                    target_coord,
                    tag="interconnect_lines",
                    geometry=geom,
                    length=geom.length,
                )
    return G


def connect_network_components(
    G, interconnect_dist=2.0, edge_tags=["theme_lines", "ref_lines"]
):
    components = list(nx.connected_components(G))
    if len(components) <= 1:
        return G

    allowed_nodes = set()
    for u, v, d in G.edges(data=True):
        if d.get("tag") in edge_tags:
            allowed_nodes.add(u)
            allowed_nodes.add(v)

    component_data = []
    for i, comp in enumerate(components):
        ref_in_comp = [n for n in comp if n in allowed_nodes]
        if ref_in_comp:
            comp_points = [Point(n) for n in ref_in_comp]
            component_data.append(
                {
                    "id": i,
                    "nodes": ref_in_comp,
                    "points": comp_points,
                    "tree": STRtree(comp_points),
                }
            )

    for i in range(len(component_data)):
        comp_a = component_data[i]
        if not comp_a["points"]:
            continue
        best_dist = float("inf")
        best_connection = None

        for j in range(len(component_data)):
            if i == j:
                continue
            comp_b = component_data[j]
            for idx_a, pt_a in enumerate(comp_a["points"]):
                idx_b = comp_b["tree"].nearest(pt_a)
                pt_b = comp_b["points"][idx_b]
                dist = pt_a.distance(pt_b)
                if dist < best_dist and dist <= interconnect_dist:
                    best_dist = dist
                    best_connection = (comp_a["nodes"][idx_a], comp_b["nodes"][idx_b])
                    if best_dist == 0:
                        break
            if best_dist == 0:
                break

        if best_connection:
            u, v = best_connection
            if not G.has_edge(u, v):
                geom = LineString([u, v])
                G.add_edge(
                    u,
                    v,
                    tag="component_interconnect",
                    geometry=geom,
                    length=geom.length,
                )
                logging.debug(
                    f"Component {i} connected to another component (distance: {best_dist:.3f})"
                )
    return G


def connect_unconnected_greedy(G, max_spatial_dist=50, detour_ratio=3.0):
    dead_ends = [n for n, d in G.degree() if d == 1]
    target_edges = []
    for u, v, data in G.edges(data=True):
        if "geometry" in data:
            if u in dead_ends or v in dead_ends:
                target_edges.append((u, v, data["geometry"]))

    all_candidates = []
    for node_id in dead_ends:
        node_point = Point(node_id)
        for u, v, edge_geom in target_edges:
            if node_id == u or node_id == v:
                continue
            p1, p2 = nearest_points(node_point, edge_geom)
            spatial_dist = p1.distance(p2)
            if spatial_dist < max_spatial_dist:
                all_candidates.append(
                    {
                        "from_node": node_id,
                        "to_edge": (u, v),
                        "projection_point": p2,
                        "spatial_dist": spatial_dist,
                    }
                )

    all_candidates.sort(key=lambda x: x["spatial_dist"])
    for link in all_candidates:
        u, v = link["to_edge"]
        dead_end_node = link["from_node"]
        proj_p = link["projection_point"]
        spatial_dist = link["spatial_dist"]
        try:
            current_network_dist = nx.shortest_path_length(
                G, source=dead_end_node, target=u, weight="length"
            )
        except nx.NetworkXNoPath:
            current_network_dist = float("inf")

        if current_network_dist / spatial_dist > detour_ratio:
            pseudo_node_id = proj_p.coords[0]
            G.add_node(pseudo_node_id, tag="pseudo_connect_vertex")
            if G.has_edge(u, v):
                G.remove_edge(u, v)
                line_up = LineString([u, pseudo_node_id])
                G.add_edge(u, pseudo_node_id, geometry=line_up, length=line_up.length)
                line_pv = LineString([pseudo_node_id, v])
                G.add_edge(pseudo_node_id, v, geometry=line_pv, length=line_pv.length)
            else:
                continue
            line_connect = LineString([dead_end_node, pseudo_node_id])
            G.add_edge(
                dead_end_node,
                pseudo_node_id,
                geometry=line_connect,
                length=line_connect.length,
                tag="connect_unconnected",
            )
    return G


def connect_ref_points_to_nearest(G, ref_points, k_neighbors=2):
    if not ref_points or ref_points.is_empty:
        return G
    tag = "ref_points"
    for pt in ref_points.geoms:
        p_coord = pt.coords[0]
        if p_coord not in G:
            G.add_node(p_coord, tag=tag)
        else:
            G.nodes[p_coord].update({"tag": tag})

    ref_nodes = [n for n, d in G.nodes(data=True) if d.get("tag") == tag]
    if len(ref_nodes) <= 1:
        return G

    coords = np.array(ref_nodes, dtype=float)
    n = len(ref_nodes)
    k = min(k_neighbors, len(ref_nodes) - 1)
    if k <= 0:
        return G

    max_matrix_nodes = 4000
    if n <= max_matrix_nodes:
        sq_norms = np.einsum("ij,ij->i", coords, coords)
        d2_matrix = sq_norms[:, None] + sq_norms[None, :] - 2.0 * coords.dot(coords.T)
        np.fill_diagonal(d2_matrix, np.inf)

        for i, u in enumerate(ref_nodes):
            d2 = d2_matrix[i]
            nearest_idx = np.argpartition(d2, k)[:k]
            nearest_idx = nearest_idx[np.argsort(d2[nearest_idx])]
            for j in nearest_idx:
                v = ref_nodes[int(j)]
                if not G.has_edge(u, v):
                    geom = LineString([Point(u), Point(v)])
                    G.add_edge(
                        u,
                        v,
                        tag="ref_points_connection",
                        geometry=geom,
                        length=geom.length,
                    )
    else:
        for i, u in enumerate(ref_nodes):
            diff = coords - coords[i]
            d2 = diff[:, 0] ** 2 + diff[:, 1] ** 2
            d2[i] = np.inf
            nearest_idx = np.argpartition(d2, k)[:k]
            nearest_idx = nearest_idx[np.argsort(d2[nearest_idx])]
            for j in nearest_idx:
                v = ref_nodes[int(j)]
                if not G.has_edge(u, v):
                    geom = LineString([Point(u), Point(v)])
                    G.add_edge(
                        u,
                        v,
                        tag="ref_points_connection",
                        geometry=geom,
                        length=geom.length,
                    )
    return G


def reconstruct_graph_for_intersections(G, tolerance=1e-6):
    edge_items = [(u, v, d) for u, v, d in G.edges(data=True) if "geometry" in d]
    if not edge_items:
        return G

    edge_geoms = [d["geometry"] for u, v, d in edge_items]
    tree = STRtree(edge_geoms)
    has_intersections = False
    for i, geom in enumerate(edge_geoms):
        possible_matches = tree.query(geom, predicate="intersects")
        for match_idx in possible_matches:
            if i >= match_idx:
                continue
            if edge_geoms[i].crosses(edge_geoms[match_idx]):
                has_intersections = True
                break
        if has_intersections:
            break

    if not has_intersections:
        return G
    return _perform_reconstruction(G, edge_geoms, tolerance)


def _perform_reconstruction(G, old_edge_geoms, tolerance):
    old_nodes_data = {n: d for n, d in G.nodes(data=True)}
    old_node_ids = [n for n, d in G.nodes(data=True)]
    old_node_geoms = [Point(n) for n in old_node_ids]
    node_tree = STRtree(old_node_geoms)
    old_edges_list = [d for u, v, d in G.edges(data=True) if "geometry" in d]
    edge_tree = STRtree(old_edge_geoms)
    noded_output = unary_union(old_edge_geoms)
    new_segments = (
        noded_output.geoms if hasattr(noded_output, "geoms") else [noded_output]
    )
    new_G = nx.Graph()

    for seg in new_segments:
        segment_node_ids = []
        for coord in [seg.coords[0], seg.coords[-1]]:
            p = Point(coord)
            idx = node_tree.nearest(p)
            if idx is not None:
                nearest_geom = old_node_geoms[idx]
                if p.distance(nearest_geom) < tolerance:
                    orig_id = old_node_ids[idx]
                    if orig_id not in new_G:
                        new_G.add_node(orig_id, **old_nodes_data[orig_id])
                    segment_node_ids.append(orig_id)
                    continue
            new_id = coord
            if new_id not in new_G:
                new_G.add_node(new_id, geometry=p, tag="pseudo_intersection_vertex")
            segment_node_ids.append(new_id)

        midpoint = seg.interpolate(0.5, normalized=True)
        edge_idx = edge_tree.nearest(midpoint)
        if edge_idx is not None:
            new_edge_data = old_edges_list[edge_idx].copy()
            new_edge_data["geometry"] = seg
            new_edge_data["length"] = seg.length
            new_G.add_edge(segment_node_ids[0], segment_node_ids[1], **new_edge_data)
    return new_G


def _candidate_nodes_within_tolerance(
    point, nodes, tolerance, node_points_tree=None, node_set=None
):
    if tolerance is None or tolerance <= 0:
        return []
    node_list = list(nodes)
    if not node_list:
        return []
    if node_points_tree is None:
        node_points = [Point(n) for n in node_list]
        tree = STRtree(node_points)
    else:
        tree = node_points_tree
    if node_set is None:
        node_set = set(node_list)
    try:
        hits = tree.query(point, predicate="dwithin", distance=tolerance)
    except Exception:
        hits = []
    candidates = []
    for h in hits:
        if isinstance(h, (int, np.integer)):
            candidates.append(node_list[int(h)])
        else:
            coord = (h.x, h.y)
            if coord in node_set:
                candidates.append(coord)
    seen = set()
    unique_candidates = []
    for c in candidates:
        if c not in seen:
            seen.add(c)
            unique_candidates.append(c)
    return unique_candidates


def _is_pseudo_graph_node(graph, node_id):
    tag = str(graph.nodes[node_id].get("tag", ""))
    return tag.startswith("pseudo_")


def _select_network_node_by_snap_strategy(
    point,
    graph,
    snap_strategy=SnapStrategy.NO_PREFERENCE,
    tolerance=None,
    angle_threshold_degrees=150.0,
    endnode_tag="is_reference_line_end",
    node_list=None,
    node_points_tree=None,
):
    node_list = list(graph.nodes) if node_list is None else list(node_list)
    if not node_list:
        return None
    node_set = set(node_list)
    candidates = _candidate_nodes_within_tolerance(
        point,
        node_list,
        tolerance,
        node_points_tree=node_points_tree,
        node_set=node_set,
    )
    real_candidates = [n for n in candidates if not _is_pseudo_graph_node(graph, n)]
    real_nodes = [n for n in node_list if not _is_pseudo_graph_node(graph, n)]

    if snap_strategy == SnapStrategy.PREFER_ENDS_AND_ANGLES:
        ranking_pool = real_candidates if real_candidates else candidates
        if not ranking_pool:
            if real_nodes:
                return nearest_node(point, real_nodes)
            return nearest_node(point, node_list)

        end_nodes = {n for n in ranking_pool if graph.nodes[n].get(endnode_tag, False)}
        angle_nodes = {
            n
            for n in ranking_pool
            if n not in end_nodes
            and _is_angle_graph_node(graph, n, angle_threshold_degrees)
        }

        def _priority(node):
            if node in end_nodes:
                rank = 0
            elif node in angle_nodes:
                rank = 1
            else:
                rank = 2
            return rank, float(Point(node).distance(point))
        return min(ranking_pool, key=_priority)

    if snap_strategy == SnapStrategy.ONLY_VERTICES:
        if real_candidates:
            return min(real_candidates, key=lambda n: Point(n).distance(point))
        if real_nodes:
            return nearest_node(point, real_nodes)
        if candidates:
            return min(candidates, key=lambda n: Point(n).distance(point))
        return nearest_node(point, node_list)

    if snap_strategy == SnapStrategy.PREFER_VERTICES:
        if real_candidates:
            return min(real_candidates, key=lambda n: Point(n).distance(point))
        if candidates:
            return min(candidates, key=lambda n: Point(n).distance(point))
        if real_nodes:
            return nearest_node(point, real_nodes)
        return nearest_node(point, node_list)

    if snap_strategy == SnapStrategy.NO_PREFERENCE:
        if candidates:
            return min(candidates, key=lambda n: Point(n).distance(point))
        return nearest_node(point, node_list)

    if candidates:
        return min(candidates, key=lambda n: Point(n).distance(point))
    return nearest_node(point, node_list)


def build_custom_network(
    input_geometry,
    theme_multiline,
    reference,
    reference_intersection,
    relevant_distance,
    gap_threshold=0.1,
    snap_dist=0.01,
):
    G = nx.Graph()

    ref_multiline, ref_points = multilinestring_multipoint_from_reference_intersection(
        reference_intersection
    )
    theme_points = get_thematic_points(input_geometry, reference_intersection)

    pseudo_ref_coords = get_pseudo_coords(reference_intersection, reference)
    pseudo_theme_coords = get_pseudo_coords(theme_multiline, input_geometry)
    ref_end_coords = get_end_coords(reference, ref_multiline)
    ref_endnode_tag = "is_reference_line_end"

    srtree_theme_lines, edge_mapping_theme_lines = _multilinestring_to_edges(
        G,
        theme_multiline,
        pseudo_coords=pseudo_theme_coords,
        node_tag="theme_vertex",
        pseudonode_tag="pseudo_theme_vertex",
        edge_tag="theme_lines",
    )
    srtree_ref_lines, edge_mapping_ref_lines = _multilinestring_to_edges(
        G,
        ref_multiline,
        pseudo_coords=pseudo_ref_coords,
        node_tag="ref_vertex",
        pseudonode_tag="pseudo_ref_vertex",
        edge_tag="ref_lines",
        end_coords=ref_end_coords,
        endnode_tag=ref_endnode_tag,
    )

    G = add_pseudonodes_to_ref_line(
        G,
        edge_mapping_ref_lines,
        ref_multiline,
        relevant_distance,
        snap_dist,
        srtree_ref_lines,
        theme_points,
    )

    connect_ref_points_to_nearest(G, ref_points, 2)

    G = connect_network_gaps(
        G=G,
        gap_dist=gap_threshold,
        snap_dist=snap_dist,
        merge_nodes=False,
    )

    G = connect_theme_and_reference(
        G=G,
        gap_dist=gap_threshold,
        interconnect_dist=2 * relevant_distance,
    )

    G = connect_network_components(
        G,
        interconnect_dist=5 * relevant_distance,
        edge_tags=["ref_lines", "ref_points_connection"],
    )
    remove_pseudonodes(G, tag="pseudo_ref_vertex", target_degrees=[1])
    remove_pseudonodes(G, tag="pseudo_theme_vertex", target_degrees=[2])
    max_unconnected_distance = 50
    unconnected_distance = 10 * relevant_distance
    if unconnected_distance > max_unconnected_distance:
        unconnected_distance = max_unconnected_distance
    unconnected_distance = 50
    G = connect_unconnected_greedy(
        G, max_spatial_dist=unconnected_distance, detour_ratio=5
    )
    remove_pseudonodes(G, tag="pseudo_theme_vertex", target_degrees=[1, 2])
    G = reconstruct_graph_for_intersections(G)
    return G


def find_best_circle_path(graph, geom_to_process, max_total_combis=1000):
    min_dist = inf
    best_cycle_line = None
    edges = [data["geometry"] for u, v, data in graph.edges(data=True)]
    polygonized = polygonize(edges)
    if isinstance(polygonized, GeometryCollection):
        individual_polygons = list(polygonized.geoms)
    else:
        individual_polygons = list(polygonized)

    if not individual_polygons:
        return None
    if len(individual_polygons) == 1:
        poly = individual_polygons[0]
        return poly.exterior if hasattr(poly, "exterior") else None

    num_to_include = 0
    for n in range(1, len(individual_polygons) + 1):
        if (2**n) - 1 > max_total_combis:
            break
        num_to_include = n
    if len(individual_polygons) > num_to_include:
        combined = safe_unary_union(individual_polygons)
        if isinstance(combined, Polygon):
            return combined.exterior if hasattr(combined, "exterior") else None
        elif isinstance(combined, (MultiPolygon, GeometryCollection)):
            largest_poly = max(combined.geoms, key=lambda p: p.area)
            return largest_poly.exterior if hasattr(largest_poly, "exterior") else None
        else:
            return None

    for i in range(1, len(individual_polygons) + 1):
        for combi in combinations(individual_polygons, i):
            combined_poly = safe_unary_union(combi)
            if hasattr(combined_poly, "exterior"):
                cycle_line = combined_poly.exterior
                vertex_distance = total_vertex_distance(cycle_line, geom_to_process)
                if vertex_distance < min_dist:
                    min_dist = vertex_distance
                    best_cycle_line = cycle_line
                    if min_dist < 1e-7:
                        return best_cycle_line
    return best_cycle_line


def longest_linestring_from_multilinestring(multilinestring):
    if multilinestring is None or multilinestring.is_empty:
        return GeometryCollection()
    if isinstance(multilinestring, LineString):
        return multilinestring
    if not isinstance(multilinestring, MultiLineString):
        logging.warning(
            "Multilinestring expected. Other type detected; empty geometry returned"
        )
        return GeometryCollection()

    graph = Graph()
    _multilinestring_to_edges(
        graph,
        multilinestring,
        node_tag="",
        edge_tag="",
        pseudo_coords=[],
        pseudonode_tag="",
    )

    if nx.is_forest(graph):
        best_path = []
        best_length = 0.0
        for component_nodes in nx.connected_components(graph):
            if not component_nodes or len(component_nodes) == 1:
                continue
            sub = graph.subgraph(component_nodes)
            start = next(iter(component_nodes))
            dist_start = nx.single_source_dijkstra_path_length(
                sub, start, weight="length"
            )
            far_a = max(dist_start, key=dist_start.get)
            dist_a = nx.single_source_dijkstra_path_length(sub, far_a, weight="length")
            far_b = max(dist_a, key=dist_a.get)
            path = nx.shortest_path(sub, source=far_a, target=far_b, weight="length")
            length = float(dist_a.get(far_b, 0.0))
            if length > best_length:
                best_length = length
                best_path = path
        if len(best_path) >= 2:
            return LineString(best_path)

    longest_path = []
    max_length = 0
    for source in graph.nodes:
        for target in graph.nodes:
            if source != target:
                try:
                    path = nx.shortest_path(
                        graph, source=source, target=target, weight="weight"
                    )
                    length = sum(
                        LineString([path[i], path[i + 1]]).length
                        for i in range(len(path) - 1)
                    )
                    if length > max_length:
                        max_length = length
                        longest_path = path
                except nx.NetworkXNoPath:
                    continue
    return LineString(longest_path)


def find_best_path_in_network(
    geom_to_process,
    graph,
    cutoff=1000,
    snap_strategy=SnapStrategy.NO_PREFERENCE,
    tolerance=None,
    angle_threshold_degrees=150.0,
):
    start_point = Point(geom_to_process.coords[0])
    end_point = Point(geom_to_process.coords[-1])
    node_list = list(graph.nodes)
    node_points_tree = None
    if tolerance is not None and tolerance > 0 and node_list:
        node_points_tree = STRtree([Point(n) for n in node_list])

    start_node = _select_network_node_by_snap_strategy(
        start_point,
        graph,
        snap_strategy=snap_strategy,
        tolerance=tolerance,
        angle_threshold_degrees=angle_threshold_degrees,
        node_list=node_list,
        node_points_tree=node_points_tree,
    )
    end_node = _select_network_node_by_snap_strategy(
        end_point,
        graph,
        snap_strategy=snap_strategy,
        tolerance=tolerance,
        angle_threshold_degrees=angle_threshold_degrees,
        node_list=node_list,
        node_points_tree=node_points_tree,
    )

    if start_node is None or end_node is None:
        return None
    if start_node == end_node:
        return find_best_circle_path(graph, geom_to_process)

    path_found = nx.has_path(graph, start_node, end_node)
    logging.debug("Path detected? - " + str(path_found))
    if not path_found:
        graph = connect_network_components(
            graph,
            50,
            edge_tags=["theme_lines", "ref_lines", "interconnect", "gap_closure"],
        )

    all_paths_generator = nx.all_simple_paths(
        graph, source=start_node, target=end_node, cutoff=cutoff
    )

    min_dist = inf
    best_line = None
    for i, path in enumerate(all_paths_generator):
        if i > cutoff:
            logging.warning(f"max paths tested while searching for geometry: {cutoff}")
            break
        if len(path) < 2:
            continue
        try:
            line = LineString(path)
            dist = total_vertex_distance(line, geom_to_process)
            if dist < min_dist:
                min_dist = dist
                best_line = line
                if dist == 0.0:
                    return best_line
        except Exception:
            pass
    return best_line


def bridge_with_straight_line(G, n):
    neighbors = list(G.neighbors(n))

    if len(neighbors) == 2:
        u, v = neighbors
        old_data_un = G.get_edge_data(u, n) or {}
        old_data_nv = G.get_edge_data(n, v) or {}
        geom_un = old_data_un.get("geometry")
        geom_nv = old_data_nv.get("geometry")

        if isinstance(geom_un, LineString) and isinstance(geom_nv, LineString):
            coords_un = list(geom_un.coords)
            coords_nv = list(geom_nv.coords)
            if coords_un and coords_un[-1] != n and coords_un[0] == n:
                coords_un = list(reversed(coords_un))
            if coords_nv and coords_nv[0] != n and coords_nv[-1] == n:
                coords_nv = list(reversed(coords_nv))
            merged_coords = coords_un
            if coords_nv:
                if merged_coords and coords_nv[0] == merged_coords[-1]:
                    merged_coords = merged_coords + coords_nv[1:]
                else:
                    merged_coords = merged_coords + coords_nv
            deduped = []
            for c in merged_coords:
                if not deduped or c != deduped[-1]:
                    deduped.append(c)
            if len(deduped) >= 2:
                new_geom = LineString(deduped)
            else:
                new_geom = LineString([u, v])
        else:
            new_geom = LineString([u, v])

        new_data = old_data_un.copy()
        new_data["geometry"] = new_geom
        new_data["length"] = new_geom.length
        G.add_edge(u, v, **new_data)
        G.remove_node(n)

    elif len(neighbors) < 2:
        G.remove_node(n)


def clean_pseudo_nodes_by_snap_strategy(
    G,
    snap_strategy=SnapStrategy.NO_PREFERENCE,
    distance_threshold=5.0,
    angle_threshold_degrees=150.0,
    endnode_tag="is_reference_line_end",
):
    if snap_strategy == SnapStrategy.NO_PREFERENCE:
        return G
    working_G = G.copy()

    def _build_tree(nodes):
        geoms = [Point(n) for n in nodes]
        return STRtree(geoms) if geoms else None

    def _current_real_nodes():
        return [
            n
            for n, d in working_G.nodes(data=True)
            if not str(d.get("tag", "")).startswith("pseudo_")
        ]

    def _current_pseudo_nodes():
        return [
            n
            for n, d in working_G.nodes(data=True)
            if str(d.get("tag", "")).startswith("pseudo_")
        ]

    def _is_angle_node(node_id):
        neighbors = list(working_G.neighbors(node_id))
        if len(neighbors) < 2:
            return False
        center = np.array([float(node_id[0]), float(node_id[1])], dtype=float)
        vectors = []
        for nb in neighbors:
            vec = np.array([float(nb[0]), float(nb[1])], dtype=float) - center
            norm = np.linalg.norm(vec)
            if norm > 0:
                vectors.append(vec / norm)
        if len(vectors) < 2:
            return False
        min_angle = 180.0
        for i in range(len(vectors)):
            for j in range(i + 1, len(vectors)):
                angle = _angle_between_vectors_degrees(vectors[i], vectors[j])
                if angle < min_angle:
                    min_angle = angle
        return min_angle <= angle_threshold_degrees

    def _remove_nodes_near(target_nodes, candidate_nodes, protected_nodes=None):
        tree = _build_tree(target_nodes)
        if tree is None:
            return
        protected_nodes = protected_nodes or set()
        for n in candidate_nodes:
            if n not in working_G or n in protected_nodes:
                continue
            p_geom = Point(n)
            hits = tree.query(p_geom, predicate="dwithin", distance=distance_threshold)
            if len(hits) > 0:
                bridge_with_straight_line(working_G, n)

    def _remove_pseudo_nodes_near(real_nodes_subset):
        _remove_nodes_near(real_nodes_subset, _current_pseudo_nodes())

    if snap_strategy == SnapStrategy.ONLY_VERTICES:
        for n in _current_pseudo_nodes():
            if n in working_G:
                bridge_with_straight_line(working_G, n)
        return working_G

    if snap_strategy == SnapStrategy.PREFER_VERTICES:
        _remove_pseudo_nodes_near(_current_real_nodes())
        return working_G

    if snap_strategy == SnapStrategy.PREFER_ENDS_AND_ANGLES:
        real_nodes = _current_real_nodes()
        end_nodes = [n for n in real_nodes if working_G.nodes[n].get(endnode_tag, False)]
        remove_candidates = [
            n for n in real_nodes if n not in end_nodes and working_G.degree(n) == 2
        ]
        _remove_nodes_near(
            end_nodes,
            remove_candidates,
            protected_nodes=set(end_nodes),
        )

        real_nodes = _current_real_nodes()
        angle_nodes = [n for n in real_nodes if _is_angle_node(n)]
        protected_nodes = set(end_nodes) | set(angle_nodes)
        remove_candidates = [
            n
            for n in real_nodes
            if n not in protected_nodes and working_G.degree(n) == 2
        ]
        _remove_nodes_near(
            angle_nodes,
            remove_candidates,
            protected_nodes=protected_nodes,
        )

        _remove_pseudo_nodes_near(_current_real_nodes())

    return working_G

def get_non_pseudo_coords(geom1, geom2):
    coords1 = get_coordinates(geom1)
    coords2 = get_coordinates(geom2)
    # (Tuples are 'hashable', NumPy arrays not)
    set1 = set(map(tuple, coords1))
    set2 = set(map(tuple, coords2))
    pseudo = set1 - set2
    return set1 - pseudo

# def connect_unconnected2(G, max_geo_dist=50, multi_factor=3):
#     """
#     Connects end-nodes to the nearest logical node (regardless of degree)
#     provided there is a significant topological detour in the graph.
#     """
#     # 1. Identify end-nodes (starting points for the search)
#     endnodes = [n for n, d in G.degree() if d == 1]
#
#     # 2. Index ALL nodes in the graph as potential targets
#     all_node_ids = list(G.nodes())
#     all_geoms = [Point(n) for n in all_node_ids]
#     tree = STRtree(all_geoms)
#
#     edges_to_add = []
#
#     for u_node in endnodes:
#         u_geom = Point(u_node)
#
#         # 3. Find all nodes within the geographical radius
#         # Returns indices corresponding to 'all_node_ids'
#         indices = tree.query(u_geom, predicate="dwithin", distance=max_geo_dist)
#
#         best_target = None
#         min_dist = float("inf")
#
#         for idx in indices:
#             v_node = all_node_ids[idx]
#
#             if u_node == v_node:
#                 continue
#             spatial_dist = u_geom.distance(Point(v_node))
#             # 4. Check the topological 'detour'
#             try:
#                 graph_dist = nx.shortest_path_length(
#                     G, source=u_node, target=v_node, weight="length"
#                 )
#             except nx.NetworkXNoPath:
#                 # If they are in different components, the detour is infinite
#                 graph_dist = float("inf")
#
#             # Only consider connecting if it represents a significant shortcut
#             if graph_dist >= multi_factor * spatial_dist:
#                 if spatial_dist < min_dist:
#                     min_dist = spatial_dist
#                     best_target = v_node
#
#         # 5. Store the best shortcut for this specific end-node
#         if best_target:
#             line = LineString([u_geom, Point(best_target)])
#             edges_to_add.append((u_node, best_target, min_dist, line))
#
#     # 6. Add the shortcuts to the graph
#     for u, v, d, line in edges_to_add:
#         # Check to ensure we haven't already created a better connection
#         if not G.has_edge(u, v):
#             G.add_edge(
#                 u, v, geometry=line, tag="connect_unconnected", length=line.length
#             )
#             # print(f"Shortcut created: Endnode {u} -> Node {v} (degree {G.degree(v)})")
#
#     return G

# def connect_unconnected(G, max_spatial_dist=50, detour_ratio=3.0):
#     """
#     Identifies dead-end nodes that are spatially close to an edge but
#     topologically far away (a large detour).
#     """
#     dead_ends = [n for n, d in G.degree() if d == 1]
#
#     # Pre-collect edge geometries to avoid repeated lookups
#     # We store these as 'target_edges'
#     target_edges = []
#     for u, v, data in G.edges(data=True):
#         if "geometry" in data:
#             if u in dead_ends or v in dead_ends:
#                 target_edges.append((u, v, data["geometry"]))
#
#     links_to_add = []
#
#     for node_id in dead_ends:
#         node_point = Point(node_id)
#
#         best_dist = float("inf")
#         best_link = None
#
#         for u, v, edge_geom in target_edges:
#             # Skip edges that are already connected to this dead-end
#             if node_id == u or node_id == v:
#                 continue
#
#             # Find the closest point on the edge to the dead-end node
#             p1, p2 = nearest_points(node_point, edge_geom)
#             spatial_dist = p1.distance(p2)
#
#             if spatial_dist < max_spatial_dist and spatial_dist < best_dist:
#                 try:
#                     # Calculate the current detour (network distance)
#                     # We use 'u' as a proxy for the network distance to that edge
#                     network_dist = nx.shortest_path_length(
#                         G, source=node_id, target=u, weight="length"
#                     )
#
#                     if network_dist / spatial_dist > detour_ratio:
#                         best_dist = spatial_dist
#                         best_link = {
#                             "from_node": node_id,
#                             "to_edge": (u, v),
#                             "projection_point": p2,
#                             "spatial_dist": spatial_dist,
#                             "network_dist": network_dist,
#                         }
#                 except nx.NetworkXNoPath:
#                     # If no path exists at all, this connection is high priority
#                     best_link = {
#                         "from_node": node_id,
#                         "to_edge": (u, v),
#                         "projection_point": p2,
#                         "spatial_dist": spatial_dist,
#                         "network_dist": float("inf"),
#                     }
#
#         if best_link:
#             links_to_add.append(best_link)
#
#     for link in links_to_add:
#         u, v = link["to_edge"]
#         dead_end_node = link["from_node"]
#         proj_p = link["projection_point"]
#
#         # 1. Create unique ID for the pseudo node
#         pseudo_node_id = proj_p.coords[0]
#
#         # Add the node with the requested tag
#         G.add_node(pseudo_node_id, tag="pseudo_connect_vertex")
#
#         # 3. Remove original edge and split it into two
#         if G.has_edge(u, v):
#             G.remove_edge(u, v)
#
#         # Segment u -> Pseudo
#         line_up = LineString([u, pseudo_node_id])
#         G.add_edge(u, pseudo_node_id, geometry=line_up, length=line_up.length)
#
#         # Segment Pseudo -> v
#         line_pv = LineString([pseudo_node_id, v])
#         G.add_edge(pseudo_node_id, v, geometry=line_pv, length=line_pv.length)
#
#         # 4. Create the final interconnection: Dead-end -> Pseudo
#         line_connect = LineString([dead_end_node, pseudo_node_id])
#         G.add_edge(
#             dead_end_node,
#             pseudo_node_id,
#             geometry=line_connect,
#             length=line_connect.length,
#             tag="connect_unconnected",
#         )
#
#     return G

__all__ = [
    "_multilinestring_to_edges",
    "add_pseudonodes_to_ref_line",
    "connect_network_gaps",
    "connect_network_components",
    "connect_ref_points_to_nearest",
    "connect_theme_and_reference",
    "connect_unconnected_greedy",
    "get_non_pseudo_coords",
    "get_end_coords",
    "get_pseudo_coords",
    "get_thematic_points",
    "nearest_node",
    "reconstruct_graph_for_intersections",
    "remove_pseudonodes",
    "build_custom_network",
    "find_best_circle_path",
    "longest_linestring_from_multilinestring",
    "total_vertex_distance",
    "find_best_path_in_network",
    "bridge_with_straight_line",
    "clean_pseudo_nodes_by_snap_strategy",
]
