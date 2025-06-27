import networkx as nx
from shapely.geometry import LineString, Point
from shapely.ops import nearest_points
import matplotlib.pyplot as plt
import itertools

# Voorbeeldsegmenten die een onderbroken cyclus vormen
segments = [
    LineString([(0, 0), (1, 0)]),
    LineString([(1, 0), (1, 1)]),
    LineString([(1, 1), (0, 1)]),
    # ontbrekend segment: (0,1) -> (0,0)
]

# Maak een graaf van de segmenten
G = nx.Graph()
for seg in segments:
    coords = list(seg.coords)
    G.add_edge(tuple(coords[0]), tuple(coords[1]), geometry=seg)

# Zoek de verbonden componenten
components = list(nx.connected_components(G))

# Zoek alle eindpunten van de componenten
endpoints = []
for comp in components:
    for node in comp:
        if G.degree[node] == 1:
            endpoints.append(node)

# Zoek het kortste segment tussen eindpunten van verschillende componenten
min_dist = float('inf')
best_pair = None
for a, b in itertools.combinations(endpoints, 2):
    if not nx.has_path(G, a, b):
        dist = Point(a).distance(Point(b))
        if dist < min_dist:
            min_dist = dist
            best_pair = (a, b)

# Voeg het kortste verbindingssegment toe
if best_pair:
    new_seg = LineString([best_pair[0], best_pair[1]])
    G.add_edge(best_pair[0], best_pair[1], geometry=new_seg)
    segments.append(new_seg)

# Visualisatie
fig, ax = plt.subplots()
for seg in segments:
    x, y = seg.xy
    ax.plot(x, y, 'b')

# Teken knopen
for node in G.nodes:
    ax.plot(*node, 'ro')

# Teken het toegevoegde segment in rood
if best_pair:
    x, y = new_seg.xy
    ax.plot(x, y, 'r--', linewidth=2, label='Toegevoegd segment')

ax.set_aspect('equal')
ax.legend()
plt.title("Graaf met automatisch toegevoegd segment om cyclus te sluiten")
plt.show()

