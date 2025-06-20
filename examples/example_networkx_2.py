import networkx as nx
import matplotlib.pyplot as plt
from shapely.geometry import Point

# Stel een netwerk op met knopen en gewogen verbindingen
G = nx.Graph()

# Voeg knopen toe (bijvoorbeeld coördinaten)
nodes = [
    (0, 0),
    (1, 2),
    (2, 1),
    (3, 3),
    (4, 0),
    (5, 2)
]

# Voeg knopen toe aan de graph
for node in nodes:
    G.add_node(node)

# Voeg verbindingen toe met gewichten (afstand tussen punten)
edges = [
    ((0, 0), (1, 2)),
    ((0, 0), (2, 1)),
    ((1, 2), (3, 3)),
    ((2, 1), (3, 3)),
    ((2, 1), (4, 0)),
    ((3, 3), (5, 2)),
    ((4, 0), (5, 2))
]

# Voeg gewichten toe op basis van euclidische afstand
for u, v in edges:
    dist = Point(u).distance(Point(v))
    G.add_edge(u, v, weight=dist)

# Stel een reeks punten in die we in volgorde willen verbinden
route_points = [(0, 0), (3, 3), (5, 2)]

# Bepaal het beste pad dat deze punten in volgorde verbindt
full_path = []
for i in range(len(route_points) - 1):
    start = route_points[i]
    end = route_points[i + 1]
    path_segment = nx.shortest_path(G, source=start, target=end, weight='weight')
    if full_path:
        # Vermijd duplicatie van knooppunten
        full_path.extend(path_segment[1:])
    else:
        full_path.extend(path_segment)

# Visualiseer het netwerk en het pad
pos = {node: node for node in G.nodes()}
nx.draw(G, pos, with_labels=True, node_color='lightgray', edge_color='gray')
path_edges = list(zip(full_path[:-1], full_path[1:]))
nx.draw_networkx_edges(G, pos, edgelist=path_edges, edge_color='red', width=2)
nx.draw_networkx_nodes(G, pos, nodelist=full_path, node_color='red')
plt.title("Shortest Path Connecting Given Points")
plt.show()

# Print het pad
print("Beste pad dat de punten verbindt:", full_path)

