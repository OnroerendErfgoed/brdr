import networkx as nx
from shapely.geometry import LineString

from brdr.geometry_utils import bridge_with_straight_line


def test_bridge_with_straight_line_endpoint_removal_no_new_edge():
    G = nx.Graph()
    G.add_node((0.0, 0.0), tag="ref_vertex")
    G.add_node((1.0, 0.0), tag="pseudo_ref_vertex")
    G.add_edge(
        (0.0, 0.0),
        (1.0, 0.0),
        tag="ref_lines",
        geometry=LineString([(0.0, 0.0), (1.0, 0.0)]),
        length=1.0,
    )

    bridge_with_straight_line(G, (1.0, 0.0))

    assert (1.0, 0.0) not in G.nodes
    assert len(G.edges()) == 0
