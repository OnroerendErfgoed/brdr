import geopandas as gpd
import matplotlib.pyplot as plt
import networkx as nx
import osmnx as ox
from shapely.geometry import Point

# 1. Download het straatnetwerk
place_name = "Amsterdam, Netherlands"
G = ox.graph_from_place(place_name, network_type="drive")

# 2. GPS-coördinaten (lat, lon)
gps_points = [
    (52.370216, 4.895168),  # Centraal Station
    (52.3676, 4.9041),  # Rembrandtplein
    (52.3600, 4.8852),  # Museumplein
]

# 3. Vind dichtstbijzijnde knooppunten
nodes = ox.distance.nearest_nodes(
    G, X=[lon for lat, lon in gps_points], Y=[lat for lat, lon in gps_points]
)

# 4. Bereken kortste pad tussen de knooppunten
route = []
for i in range(len(nodes) - 1):
    path = nx.shortest_path(G, nodes[i], nodes[i + 1], weight="length")
    route.extend(path if i == 0 else path[1:])  # vermijd dubbele knooppunten

# 5. Visualiseer het netwerk en de route
fig, ax = ox.plot_graph_route(
    G, route, route_linewidth=4, node_size=0, bgcolor="white", show=False, close=False
)

# 6. Voeg originele GPS-punten toe
gps_geom = [Point(lon, lat) for lat, lon in gps_points]
gdf_gps = gpd.GeoDataFrame(geometry=gps_geom, crs="EPSG:4326").to_crs(G.graph["crs"])
gdf_gps.plot(ax=ax, color="red", markersize=50, label="GPS-punten")

plt.legend()
plt.title("Map-matching van GPS-traject op straatnetwerk")
plt.tight_layout()
plt.show()
