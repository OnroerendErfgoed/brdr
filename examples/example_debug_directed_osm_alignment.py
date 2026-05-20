from collections import Counter

from shapely import wkt
from shapely.affinity import translate
from shapely.geometry import LineString, MultiLineString

from brdr.aligner import Aligner
from brdr.configs import ProcessorConfig
from brdr.enums import SnapStrategy
from brdr.loader import DictLoader
from brdr.osm.loader import OSMLoader
from brdr.processor import NetworkGeometryProcessor

CRS = "EPSG:3812"
INPUT_ID = "theme_line"
RELEVANT_DISTANCE = 20
CORRECTION_DISTANCE = 0.01

THEMATIC_INPUT = wkt.loads(
    "LINESTRING (655355.8635293404 713357.4894204946, "
    "655382.1514974725 713370.1638195604, "
    "655422.1900856637 713388.3282459967, "
    "655488.0750510992 713418.2084318392, "
    "655496.4051359073 713421.809954524, "
    "655547.3701347096 713443.6012302733, "
    "655598.3437011356 713464.4391079134)"
)


def _print_header(title: str):
    print("")
    print("=" * 88)
    print(title)
    print("=" * 88)


def _line_direction_info(line: LineString, label: str):
    start = tuple(line.coords[0])
    end = tuple(line.coords[-1])
    print(f"{label} start: {start}")
    print(f"{label} end  : {end}")
    print(f"{label} length: {round(line.length, 3)} m")


def _load_osm_reference(seed_line: LineString):
    _print_header("STEP 1 - Load OSM Reference Data")
    aligner_seed = Aligner(crs=CRS)
    aligner_seed.load_thematic_data(DictLoader({INPUT_ID: seed_line}))
    loader = OSMLoader(
        osm_tags={"highway": True},
        aligner=aligner_seed,
        included_attributes=("osm_id", "element", "highway", "name", "oneway"),
        include_directional_attributes=True,
    )
    aligner_seed.load_reference_data(loader)

    features = aligner_seed.reference_data.features
    print(f"Reference feature count: {len(features)}")
    oneway_counter = Counter()
    highway_counter = Counter()
    missing_oneway = 0
    for feat in features.values():
        props = feat.properties or {}
        highway = str(props.get("highway", "")).strip().lower()
        if highway:
            highway_counter[highway] += 1
        oneway = str(props.get("oneway", "")).strip().lower()
        if oneway:
            oneway_counter[oneway] += 1
        else:
            missing_oneway += 1

    print("Top highway values:")
    for key, value in highway_counter.most_common(10):
        print(f"  - {key}: {value}")
    print("Oneway distribution:")
    for key, value in oneway_counter.items():
        print(f"  - {key}: {value}")
    print(f"Missing/empty oneway: {missing_oneway}")

    reference_dict = {key: feat.geometry for key, feat in features.items()}
    reference_props = {
        key: dict(feat.properties or {}) for key, feat in features.items()
    }
    return reference_dict, reference_props


def _run_alignment(
    reference_dict, reference_props, input_line: LineString, directed: bool
):
    mode = "DIRECTED" if directed else "UNDIRECTED"
    _print_header(f"STEP 2 - Run Alignment ({mode})")
    _line_direction_info(input_line, f"Input line ({mode})")

    config = ProcessorConfig(
        snap_strategy=SnapStrategy.NO_PREFERENCE,
        network_use_directed_graph=directed,
        network_oneway_field="oneway",
        network_oneway_forward_values=("yes", "1", "true", "forward"),
        network_oneway_reverse_values=("-1", "reverse", "backward"),
        network_allow_connector_edges_when_directed=True,
        network_directed_connector_penalty_factor=10.0,
    )
    print("Config:")
    print(f"  - network_use_directed_graph={config.network_use_directed_graph}")
    print(f"  - network_oneway_field={config.network_oneway_field}")
    print(f"  - forward_values={config.network_oneway_forward_values}")
    print(f"  - reverse_values={config.network_oneway_reverse_values}")
    print(
        "  - allow_connector_edges_when_directed="
        f"{config.network_allow_connector_edges_when_directed}"
    )
    print(
        "  - directed_connector_penalty_factor="
        f"{config.network_directed_connector_penalty_factor}"
    )

    aligner = Aligner(crs=CRS, processor=NetworkGeometryProcessor(config=config))
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
    result = out.results[INPUT_ID][RELEVANT_DISTANCE]["result"]
    if result is None:
        print("Result: None (no valid route found)")
    else:
        print(f"Result length: {round(result.length, 3)} m")
        print(f"Result WKT   : {result.wkt}")
    return result


def _compare_results(undirected_result, directed_result):
    _print_header("STEP 3 - Compare Directed vs Undirected")
    if undirected_result is None and directed_result is None:
        print("Both results are None.")
        return
    if undirected_result is None:
        print("Undirected is None, directed found a route.")
        return
    if directed_result is None:
        print("Directed is None, undirected found a route.")
        return
    symmetric_diff_len = undirected_result.symmetric_difference(directed_result).length
    overlap_len = undirected_result.intersection(directed_result).length
    print(f"Symmetric difference length: {round(symmetric_diff_len, 3)} m")
    print(f"Overlap length              : {round(overlap_len, 3)} m")
    print(
        "Observation: if directed logic works for one-way separation, "
        "you typically see less side-switching and a clear geometry difference."
    )


def _run_reverse_direction_check(reference_dict, reference_props):
    _print_header("STEP 4 - Forward vs Reverse Direction Check")
    forward_line = THEMATIC_INPUT
    reverse_line = LineString(list(THEMATIC_INPUT.coords)[::-1])
    reverse_line_shifted = translate(reverse_line, xoff=0.3, yoff=0.3)

    print("Forward input:")
    _line_direction_info(forward_line, "Forward")
    print("Reverse input (slightly shifted):")
    _line_direction_info(reverse_line_shifted, "Reverse")

    forward_directed = _run_alignment(
        reference_dict, reference_props, forward_line, directed=True
    )
    reverse_directed = _run_alignment(
        reference_dict, reference_props, reverse_line_shifted, directed=True
    )

    print("")
    print("Forward directed result is None? ", forward_directed is None)
    print("Reverse directed result is None? ", reverse_directed is None)
    print(
        "Expectation: reverse should fail more often or follow another corridor "
        "when one-way constraints are effective."
    )


def main():
    _print_header("DIRECTED OSM DEBUG - SINGLE CARRIAGEWAY ALIGNMENT")
    print("Goal: verify that directed routing stays on the correct one-way side.")
    print(f"CRS: {CRS}")
    print(f"Relevant distance: {RELEVANT_DISTANCE}")
    print(f"Correction distance: {CORRECTION_DISTANCE}")
    _line_direction_info(THEMATIC_INPUT, "Given thematic input")

    reference_dict, reference_props = _load_osm_reference(THEMATIC_INPUT)
    undirected_result = _run_alignment(
        reference_dict, reference_props, THEMATIC_INPUT, directed=False
    )
    directed_result = _run_alignment(
        reference_dict, reference_props, THEMATIC_INPUT, directed=True
    )
    _compare_results(undirected_result, directed_result)
    _run_reverse_direction_check(reference_dict, reference_props)


if __name__ == "__main__":
    main()
