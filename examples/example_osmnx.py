# Demo osmx usage
# https://medium.com/@callumjamesscoby/routing-shortest-path-analysis-in-python-and-qgis-39da55da5099
import osmnx as ox

# Choose tags (fe cycleway)
tags = {"highway": "cycleway"}
# Download cycleways in Leuven
leuven_cycleways = ox.features_from_place("Leuven, Belgium", tags)
# show rows
print(leuven_cycleways.head())
# Sava as geopackage
leuven_cycleways.to_file("output/leuven_cycleways.gpkg", driver="GPKG")


# Demo osmnx routing
G = ox.graph_from_place("Heverlee, Belgium", network_type="drive")
Gp = ox.project_graph(G, to_crs=31370)
ox.plot_graph(Gp)
ox.save_graph_geopackage(Gp, filepath="output/heverlee.gpkg")
list(Gp)

X0 = 200000, 0  # Point A
X1 = 200000, 0  # Point B
# 50,849374°  4,693125°

node = ox.nearest_nodes(Gp, X0, X1)

X0_node = [i for i, x in enumerate(list(Gp)) if x == node[0]]
X1_node = [i for i, x in enumerate(list(Gp)) if x == node[1]]

print(X0_node, X1_node)

orig = list(Gp)[X0_node[0]]
dest = list(Gp)[X1_node[0]]

shortest_distance_route = ox.shortest_path(Gp, orig, dest, weight="length")

fig, ax = ox.plot_graph_route(
    Gp, shortest_distance_route, route_color="r", route_linewidth=6, node_size=0
)

mutliple_routes = ox.k_shortest_paths(Gp, orig, dest, k=2, weight="length")

fig, ax = ox.plot_graph_routes(
    Gp, list(mutliple_routes), route_colors="y", route_linewidth=4, node_size=0
)
