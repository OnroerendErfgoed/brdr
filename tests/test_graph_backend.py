import pytest

from brdr.configs import ProcessorConfig
from brdr.graph_backend import (
    Graph,
    NoPath,
    _materialize_rustworkx_graph,
    all_simple_paths,
    astar_path,
    backend_name,
    has_path,
    is_backend_available,
    new_graph,
    shortest_path,
    shortest_path_length,
)


def test_new_graph_defaults_to_rustworkx_backend():
    graph = new_graph()

    assert backend_name(graph) == "rustworkx"
    assert not graph.is_directed()


def test_processor_config_defaults_to_rustworkx_backend():
    config = ProcessorConfig()

    assert config.network_graph_backend == "rustworkx"


def test_graph_copy_and_to_undirected_preserve_backend_preference():
    graph = new_graph(directed=True)
    graph.add_edge("a", "b", length=1.0)

    copied = graph.copy()
    undirected = graph.to_undirected()

    assert backend_name(copied) == "rustworkx"
    assert copied.is_directed()
    assert backend_name(undirected) == "rustworkx"
    assert not undirected.is_directed()


def test_graph_constructor_accepts_explicit_backend():
    graph = Graph(backend="networkx")

    assert backend_name(graph) == "networkx"


def test_unavailable_optional_backend_fails_early():
    backend = "igraph"
    assert not is_backend_available(backend)
    with pytest.raises(ImportError):
        new_graph(backend=backend)


def test_rustworkx_backend_uses_native_query_adapter_when_available():
    if not is_backend_available("rustworkx"):
        pytest.skip("rustworkx is not installed")

    graph = new_graph(backend="rustworkx")
    graph.add_edge("A", "B", length=2.0)
    graph.add_edge("B", "C", length=3.0)

    assert backend_name(graph) == "rustworkx"
    assert has_path(graph, "A", "C")
    assert shortest_path(graph, "A", "C", weight="length") == ["A", "B", "C"]
    assert shortest_path_length(graph, "A", "C", weight="length") == 5.0
    assert list(all_simple_paths(graph, "A", "C", cutoff=5)) == [["A", "B", "C"]]
    assert astar_path(
        graph,
        "A",
        "C",
        heuristic=lambda _a, _b: 0.0,
        weight=lambda u, v, data: data["length"],
    ) == ["A", "B", "C"]


def test_rustworkx_backend_caches_converted_graph_between_queries(monkeypatch):
    if not is_backend_available("rustworkx"):
        pytest.skip("rustworkx is not installed")

    graph = new_graph(backend="rustworkx")
    graph.add_edge("A", "B", length=2.0)
    graph.add_edge("B", "C", length=3.0)

    calls = {"count": 0}

    def counted_materialize(base_graph):
        calls["count"] += 1
        return _materialize_rustworkx_graph(base_graph)

    monkeypatch.setattr(
        "brdr.graph_backend._materialize_rustworkx_graph",
        counted_materialize,
    )

    assert shortest_path(graph, "A", "C", weight="length") == ["A", "B", "C"]
    assert shortest_path_length(graph, "A", "C", weight="length") == 5.0
    assert has_path(graph, "A", "C")
    assert calls["count"] == 1


def test_rustworkx_backend_invalidates_cached_conversion_after_graph_mutation(monkeypatch):
    if not is_backend_available("rustworkx"):
        pytest.skip("rustworkx is not installed")

    graph = new_graph(backend="rustworkx")
    graph.add_edge("A", "B", length=2.0)
    graph.add_edge("B", "C", length=3.0)

    calls = {"count": 0}

    def counted_materialize(base_graph):
        calls["count"] += 1
        return _materialize_rustworkx_graph(base_graph)

    monkeypatch.setattr(
        "brdr.graph_backend._materialize_rustworkx_graph",
        counted_materialize,
    )

    assert shortest_path(graph, "A", "C", weight="length") == ["A", "B", "C"]
    graph.add_edge("C", "D", length=1.0)
    assert shortest_path(graph, "A", "D", weight="length") == ["A", "B", "C", "D"]
    assert calls["count"] == 2


def test_rustworkx_backend_raises_nopath_for_unreachable_target():
    if not is_backend_available("rustworkx"):
        pytest.skip("rustworkx is not installed")

    graph = new_graph(backend="rustworkx", directed=True)
    graph.add_edge("A", "B", length=1.0)

    with pytest.raises(NoPath):
        shortest_path_length(graph, "B", "A", weight="length")

    with pytest.raises(NoPath):
        shortest_path(graph, "B", "A", weight="length")
