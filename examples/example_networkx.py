import matplotlib.pyplot as plt
import networkx as nx
from shapely.geometry import LineString, Point

# Voorbeeldset van LineStrings
lines = [
    LineString([(0, 0), (1, 1)]),
    LineString([(1, 1), (2, 2)]),
    LineString([(2, 2), (3, 1)]),
    LineString([(1, 1), (1, 2)]),
    LineString([(1, 2), (2, 3)])
]

# Bouw een NetworkX-graph op basis van de lijnen
G = nx.Graph()
for line in lines:
    coords = list(line.coords)
    for i in range(len(coords) - 1):
        u = coords[i]
        v = coords[i + 1]
        dist = Point(u).distance(Point(v))
        G.add_edge(u, v, weight=dist)

# Definieer start- en eindpunt
start = (0, 0)
end = (2, 3)

# Bereken het kortste pad
shortest_path = nx.shortest_path(G, source=start, target=end, weight='weight')

# Visualiseer het netwerk en het kortste pad
pos = {node: node for node in G.nodes()}
nx.draw(G, pos, with_labels=True, node_color='lightblue', edge_color='gray')
path_edges = list(zip(shortest_path[:-1], shortest_path[1:]))
nx.draw_networkx_edges(G, pos, edgelist=path_edges, edge_color='red', width=2)
plt.title("Shortest Path in Custom LineString Network")
plt.show()

# Print het pad
print("Shortest path:", shortest_path)
