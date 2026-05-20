from shapely.affinity import translate
from shapely.geometry import LineString, MultiLineString

from brdr.aligner import Aligner
from brdr.configs import ProcessorConfig
from brdr.enums import SnapStrategy
from brdr.loader import DictLoader
from brdr.osm.loader import OSMLoader
from brdr.processor import NetworkGeometryProcessor

# Example settings
CRS = "EPSG:31370"
RELEVANT_DISTANCE = 8.0
CORRECTION_DISTANCE = 0.5
INPUT_ID = "input_line"

# Seed geometry in a dense street area (Brussels center, Lambert72)
# Used to define the OSM download bbox.
SEED_THEMATIC = LineString(
    [
        (148900.0, 170650.0),
        (149350.0, 170650.0),
    ]
)


def _iter_forward_oneway_lines(reference_features):
    forward_values = {"yes", "1", "true", "forward"}
    for feat_id, feat in reference_features.items():
        props = feat.properties or {}
        oneway = str(props.get("oneway", "")).strip().lower()
        if oneway not in forward_values:
            continue
        geom = feat.geometry
        if isinstance(geom, LineString) and len(geom.coords) >= 2:
            yield feat_id, geom, props
        elif isinstance(geom, MultiLineString):
            for line in geom.geoms:
                if len(line.coords) >= 2:
                    yield feat_id, line, props


def _run_alignment(reference_dict, reference_props, input_line, directed):
    config = ProcessorConfig(
        snap_strategy=SnapStrategy.NO_PREFERENCE,
        network_use_directed_graph=directed,
        network_oneway_field="oneway",
        network_oneway_forward_values=("yes", "1", "true", "forward"),
        network_oneway_reverse_values=("-1", "reverse", "backward"),
        network_allow_connector_edges_when_directed=True,
        network_directed_connector_penalty_factor=10.0,
    )
    processor = NetworkGeometryProcessor(config=config)
    aligner = Aligner(crs=CRS, processor=processor)
    aligner.correction_distance = CORRECTION_DISTANCE
    aligner.load_thematic_data(DictLoader({INPUT_ID: input_line}))
    aligner.load_reference_data(
        DictLoader(
            data_dict=reference_dict,
            data_dict_properties=reference_props,
            is_reference=True,
        )
    )
    out = aligner.process([RELEVANT_DISTANCE], max_workers=-1)
    return out.results[INPUT_ID][RELEVANT_DISTANCE]["result"]


def main():
    # 1) Load OSM reference with attributes (incl. `oneway`)
    aligner_seed = Aligner(crs=CRS)
    aligner_seed.load_thematic_data(DictLoader({INPUT_ID: SEED_THEMATIC}))
    aligner_seed.load_reference_data(
        OSMLoader(
            osm_tags={"highway": True},
            aligner=aligner_seed,
            included_attributes=("osm_id", "element", "oneway"),
            include_directional_attributes=True,
        )
    )

    reference_dict = {
        key: feat.geometry for key, feat in aligner_seed.reference_data.features.items()
    }
    reference_props = {
        key: dict(feat.properties or {})
        for key, feat in aligner_seed.reference_data.features.items()
    }

    # 2) Find a forward-oneway line where directed and undirected differ
    selected = None
    for feat_id, oneway_line, props in _iter_forward_oneway_lines(
        aligner_seed.reference_data.features
    ):
        if oneway_line.length < 20:
            continue
        # Input runs opposite to one-way direction and slightly shifted (= not aligned)
        reverse_line = LineString(list(oneway_line.coords)[::-1])
        input_line = translate(reverse_line, xoff=0.7, yoff=0.7)

        undirected_result = _run_alignment(
            reference_dict, reference_props, input_line, directed=False
        )
        directed_result = _run_alignment(
            reference_dict, reference_props, input_line, directed=True
        )

        if undirected_result is None or directed_result is None:
            continue

        overlap_buffer = oneway_line.buffer(0.75)
        undirected_overlap = undirected_result.intersection(overlap_buffer).length
        directed_overlap = directed_result.intersection(overlap_buffer).length
        if undirected_overlap > directed_overlap + 5:
            selected = {
                "feat_id": feat_id,
                "oneway": props.get("oneway"),
                "input_line": input_line,
                "oneway_line": oneway_line,
                "undirected_result": undirected_result,
                "directed_result": directed_result,
                "undirected_overlap": undirected_overlap,
                "directed_overlap": directed_overlap,
            }
            break

    if selected is None:
        raise RuntimeError(
            "No clear directed-vs-undirected case found in this bbox. "
            "Try adjusting SEED_THEMATIC or RELEVANT_DISTANCE."
        )

    # 3) Report effect
    print("Selected OSM feature:", selected["feat_id"])
    print("oneway:", selected["oneway"])
    print("Input length:", round(selected["input_line"].length, 3))
    print("Undirected result length:", round(selected["undirected_result"].length, 3))
    print("Directed result length:", round(selected["directed_result"].length, 3))
    print(
        "Overlap with forward-oneway corridor (m):",
        "undirected=",
        round(selected["undirected_overlap"], 3),
        "directed=",
        round(selected["directed_overlap"], 3),
    )
    print("Undirected result WKT:", selected["undirected_result"].wkt)
    print("Directed result WKT:", selected["directed_result"].wkt)


if __name__ == "__main__":
    main()
