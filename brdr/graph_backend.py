"""Graph backend facade and backend-carrying graph wrapper.

The backend preference is carried explicitly through graph construction and
processor config so path-search code can default to `rustworkx` while keeping
the NetworkX-compatible graph API used across the codebase.
"""

from __future__ import annotations

from importlib.util import find_spec
from typing import Any

import networkx as nx

try:
    import rustworkx as rx
except ImportError:  # pragma: no cover - optional dependency
    rx = None

NoPath = nx.NetworkXNoPath
NodeNotFound = nx.NodeNotFound
GraphError = nx.NetworkXError

DEFAULT_BACKEND = "rustworkx"
SUPPORTED_BACKENDS = frozenset({"networkx", "rustworkx", "igraph"})
OPTIONAL_BACKEND_MODULES = {
    "rustworkx": "rustworkx",
    "igraph": "igraph",
}
_MUTATING_GRAPH_METHODS = frozenset(
    {
        "add_edge",
        "add_edges_from",
        "add_node",
        "add_nodes_from",
        "add_weighted_edges_from",
        "clear",
        "clear_edges",
        "remove_edge",
        "remove_edges_from",
        "remove_node",
        "remove_nodes_from",
        "update",
    }
)


def normalize_backend_name(backend: str | None) -> str:
    backend_name = (backend or DEFAULT_BACKEND).strip().lower()
    if backend_name not in SUPPORTED_BACKENDS:
        supported = ", ".join(sorted(SUPPORTED_BACKENDS))
        raise ValueError(
            f"Unsupported graph backend '{backend_name}'. Supported backends: {supported}."
        )
    return backend_name


def is_backend_available(backend: str | None) -> bool:
    backend_name = normalize_backend_name(backend)
    module_name = OPTIONAL_BACKEND_MODULES.get(backend_name)
    if module_name is None:
        return True
    return find_spec(module_name) is not None


def require_backend_available(backend: str | None) -> str:
    backend_name = normalize_backend_name(backend)
    if is_backend_available(backend_name):
        return backend_name
    module_name = OPTIONAL_BACKEND_MODULES[backend_name]
    raise ImportError(
        f"Graph backend '{backend_name}' requires optional dependency '{module_name}', "
        "but it is not installed in the active environment."
    )


class BackendGraph:
    """Thin wrapper around an underlying graph object.

    The wrapper preserves a backend preference across graph-producing methods
    such as ``copy()``, ``subgraph()`` and ``to_undirected()`` while delegating
    the full graph API to the wrapped object.
    """

    def __init__(
        self,
        *,
        directed: bool = False,
        backend: str = DEFAULT_BACKEND,
        graph=None,
    ):
        self.backend = require_backend_available(backend)
        self.directed = bool(directed)
        self._graph = graph if graph is not None else _new_networkx_graph(directed)
        self._rx_cached_graph = None
        self._rx_cached_node_to_index = None
        self._rx_cached_index_to_node = None

    def __getattr__(self, name: str):
        attr = getattr(self._graph, name)
        if not callable(attr) or name not in _MUTATING_GRAPH_METHODS:
            return attr

        def wrapped(*args, **kwargs):
            result = attr(*args, **kwargs)
            self._invalidate_backend_caches()
            return result

        return wrapped

    def __contains__(self, item) -> bool:
        return item in self._graph

    def __iter__(self):
        return iter(self._graph)

    def __len__(self) -> int:
        return len(self._graph)

    def __getitem__(self, key):
        return self._graph[key]

    def copy(self):
        return self.__class__(
            backend=self.backend,
            graph=self._graph.copy(),
        )

    def to_undirected(self):
        return Graph(
            backend=self.backend,
            graph=self._graph.to_undirected(),
        )

    def subgraph(self, nodes):
        return self.__class__(
            backend=self.backend,
            graph=self._graph.subgraph(nodes),
        )

    def edge_subgraph(self, edges):
        return self.__class__(
            backend=self.backend,
            graph=self._graph.edge_subgraph(edges),
        )

    def _invalidate_backend_caches(self):
        self._rx_cached_graph = None
        self._rx_cached_node_to_index = None
        self._rx_cached_index_to_node = None

    def unwrap(self):
        return self._graph


class Graph(BackendGraph):
    def __init__(self, *, backend: str = DEFAULT_BACKEND, graph=None):
        super().__init__(directed=False, backend=backend, graph=graph)


class DiGraph(BackendGraph):
    def __init__(self, *, backend: str = DEFAULT_BACKEND, graph=None):
        super().__init__(directed=True, backend=backend, graph=graph)


def _new_networkx_graph(directed: bool):
    return nx.DiGraph() if directed else nx.Graph()


def unwrap_graph(graph):
    return graph.unwrap() if isinstance(graph, BackendGraph) else graph


def backend_name(graph) -> str:
    if isinstance(graph, BackendGraph):
        return graph.backend
    return DEFAULT_BACKEND


def _use_rustworkx_backend(graph) -> bool:
    return backend_name(graph) == "rustworkx" and rx is not None


def _materialize_rustworkx_graph(base_graph):
    rx_graph = rx.PyDiGraph() if base_graph.is_directed() else rx.PyGraph()
    node_to_index = {}
    index_to_node = {}

    for node in base_graph.nodes:
        idx = rx_graph.add_node(node)
        node_to_index[node] = idx
        index_to_node[idx] = node

    for u, v, attrs in base_graph.edges(data=True):
        rx_graph.add_edge(
            node_to_index[u],
            node_to_index[v],
            {"u": u, "v": v, "attrs": attrs},
        )
    return rx_graph, node_to_index, index_to_node


def _build_rustworkx_graph(graph):
    if isinstance(graph, BackendGraph):
        if graph._rx_cached_graph is None:
            (
                graph._rx_cached_graph,
                graph._rx_cached_node_to_index,
                graph._rx_cached_index_to_node,
            ) = _materialize_rustworkx_graph(graph.unwrap())
        return (
            graph._rx_cached_graph,
            graph._rx_cached_node_to_index,
            graph._rx_cached_index_to_node,
        )
    return _materialize_rustworkx_graph(unwrap_graph(graph))


def _rustworkx_weight_fn(weight):
    if callable(weight):
        return lambda payload: weight(payload["u"], payload["v"], payload["attrs"])
    if weight is None:
        return lambda _payload: 1.0
    if isinstance(weight, str):
        return lambda payload: float(payload["attrs"].get(weight, 1.0))
    return lambda _payload: float(weight)


def _rustworkx_path_to_nodes(path, index_to_node):
    return [index_to_node[int(idx)] for idx in path]


def new_graph(*, directed: bool = False, backend: str = DEFAULT_BACKEND):
    backend_name = require_backend_available(backend)
    return DiGraph(backend=backend_name) if directed else Graph(backend=backend_name)


def connected_components(graph):
    return nx.connected_components(unwrap_graph(graph))


def weakly_connected_components(graph):
    return nx.weakly_connected_components(unwrap_graph(graph))


def is_forest(graph) -> bool:
    return nx.is_forest(unwrap_graph(graph))


def relabel_nodes(graph, mapping: dict[Any, Any], *, copy: bool = False) -> None:
    nx.relabel_nodes(unwrap_graph(graph), mapping, copy=copy)
    if isinstance(graph, BackendGraph) and not copy:
        graph._invalidate_backend_caches()


def selfloop_edges(graph):
    return nx.selfloop_edges(unwrap_graph(graph))


def has_path(graph, source, target) -> bool:
    if _use_rustworkx_backend(graph):
        rx_graph, node_to_index, _ = _build_rustworkx_graph(graph)
        if source not in node_to_index or target not in node_to_index:
            return False
        return bool(rx.has_path(rx_graph, node_to_index[source], node_to_index[target]))
    return nx.has_path(unwrap_graph(graph), source, target)


def shortest_path(graph, source, target, weight=None):
    if _use_rustworkx_backend(graph):
        rx_graph, node_to_index, index_to_node = _build_rustworkx_graph(graph)
        if source not in node_to_index or target not in node_to_index:
            raise NodeNotFound("Source or target node not found")
        try:
            path_mapping = rx.dijkstra_shortest_paths(
                rx_graph,
                node_to_index[source],
                node_to_index[target],
                weight_fn=_rustworkx_weight_fn(weight),
            )
            path = path_mapping[node_to_index[target]]
        except (KeyError, IndexError, rx.NoPathFound) as exc:
            raise NoPath(str(exc)) from exc
        return _rustworkx_path_to_nodes(path, index_to_node)
    return nx.shortest_path(
        unwrap_graph(graph), source=source, target=target, weight=weight
    )


def shortest_path_length(graph, source, target, weight=None):
    if _use_rustworkx_backend(graph):
        rx_graph, node_to_index, _ = _build_rustworkx_graph(graph)
        if source not in node_to_index or target not in node_to_index:
            raise NodeNotFound("Source or target node not found")
        lengths = rx.dijkstra_shortest_path_lengths(
            rx_graph,
            node_to_index[source],
            edge_cost_fn=_rustworkx_weight_fn(weight),
        )
        try:
            return lengths[node_to_index[target]]
        except (KeyError, IndexError) as exc:
            raise NoPath(str(exc)) from exc
    return nx.shortest_path_length(
        unwrap_graph(graph), source=source, target=target, weight=weight
    )


def single_source_dijkstra_path_length(graph, source, weight=None):
    if _use_rustworkx_backend(graph):
        rx_graph, node_to_index, index_to_node = _build_rustworkx_graph(graph)
        if source not in node_to_index:
            raise NodeNotFound("Source node not found")
        lengths = rx.dijkstra_shortest_path_lengths(
            rx_graph,
            node_to_index[source],
            edge_cost_fn=_rustworkx_weight_fn(weight),
        )
        return {index_to_node[int(idx)]: dist for idx, dist in lengths.items()}
    return nx.single_source_dijkstra_path_length(
        unwrap_graph(graph), source, weight=weight
    )


def shortest_simple_paths(graph, source, target, weight=None):
    return nx.shortest_simple_paths(
        unwrap_graph(graph), source, target, weight=weight
    )


def all_simple_paths(graph, source, target, cutoff=None):
    if _use_rustworkx_backend(graph):
        rx_graph, node_to_index, index_to_node = _build_rustworkx_graph(graph)
        if source not in node_to_index or target not in node_to_index:
            return iter(())
        paths = rx.all_simple_paths(
            rx_graph,
            node_to_index[source],
            node_to_index[target],
            cutoff=cutoff,
        )
        return (_rustworkx_path_to_nodes(path, index_to_node) for path in paths)
    return nx.all_simple_paths(unwrap_graph(graph), source, target, cutoff=cutoff)


def astar_path(graph, source, target, heuristic=None, weight=None):
    if _use_rustworkx_backend(graph):
        rx_graph, node_to_index, index_to_node = _build_rustworkx_graph(graph)
        if source not in node_to_index or target not in node_to_index:
            raise NodeNotFound("Source or target node not found")
        goal_fn = lambda node_weight: node_weight == target
        estimate_cost_fn = (
            (lambda node_weight: float(heuristic(node_weight, target)))
            if heuristic is not None
            else (lambda _node_weight: 0.0)
        )
        try:
            path = rx.astar_shortest_path(
                rx_graph,
                node_to_index[source],
                goal_fn,
                _rustworkx_weight_fn(weight),
                estimate_cost_fn,
            )
        except rx.NoPathFound as exc:
            raise NoPath(str(exc)) from exc
        return _rustworkx_path_to_nodes(path, index_to_node)
    return nx.astar_path(
        unwrap_graph(graph),
        source,
        target,
        heuristic=heuristic,
        weight=weight,
    )
