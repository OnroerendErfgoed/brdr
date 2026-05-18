import numpy as np
from shapely import from_wkt
from shapely.geometry import box

from brdr.aligner import Aligner
from brdr.be.grb.enums import GRBType
from brdr.be.grb.loader import GRBActualLoader
from brdr.configs import ProcessorConfig
from brdr.constants import EVALUATION_FIELD_NAME, PREDICTION_SCORE
from brdr.constants import RELEVANT_DISTANCE_FIELD_NAME, REMARK_FIELD_NAME
from brdr.enums import AlignerResultType, FullReferenceStrategy
from brdr.loader import DictLoader
from brdr.processor import NetworkGeometryProcessor

if __name__ == "__main__":
    """
    Example:
    - Align one polygon to GRB ADP parcels.
    - Only apply edits inside a selected processing area.
    """

    thematic_wkt = (
        "POLYGON ((173047.17180318027385511 171248.41586496974923648, "
        "173049.09147341555217281 171248.40487670945003629, "
        "173050.55786158866249025 171248.39648305167793296, "
        "173051.82306451365002431 171248.38924098529969342, "
        "173056.53338920549140312 171248.36223203572444618, "
        "173056.83458719000918791 171248.36050496983807534, "
        "173057.32397181721171364 171248.35771020868560299, "
        "173057.97769119014265016 171248.3539769696071744, "
        "173059.01052276318660006 171245.80159597989404574, "
        "173061.52002352738054469 171239.60000202123774216, "
        "173060.7997149107104633 171239.53592965667485259, "
        "173051.06421917545958422 171238.66994496964616701, "
        "173048.89839586053858511 171238.55792140372795984, "
        "173044.34415763989090919 171238.32236109697259963, "
        "173041.01759774464881048 171238.15032925445120782, "
        "173035.38723123903037049 171237.85915688565000892, "
        "173034.93621917761629447 171237.83583296000142582, "
        "173034.32331781095126644 171237.85248173048603348, "
        "173032.98882742403657176 171237.88873164856340736, "
        "173032.3350120039540343 171237.90649180306354538, "
        "173031.99542837019544095 171237.91571620738250203, "
        "173031.5773093054885976 171237.92707393763703294, "
        "173026.64902394328964874 171238.61931811177055351, "
        "173025.97929562485660426 171241.76031645329203457, "
        "173034.05359425843926147 171242.86947506395517848, "
        "173033.80399517982732505 171248.49240896999253891, "
        "173038.16392805395298637 171248.4674427515710704, "
        "173042.73404317974927835 171248.44127296973601915, "
        "173044.35234131649485789 171248.43199300387641415, "
        "173047.17180318027385511 171248.41586496974923648))"
    )

    thematic_geom = from_wkt(thematic_wkt)
    thematic = {"theme_1": thematic_geom}

    # Extent from user:
    # 173019.10, 171237.48 : 173040.45, 171248.84
    processing_area = box(173019.10, 171237.48, 173040.45, 171248.84)


    aligner = Aligner(crs="EPSG:31370")
    aligner.load_thematic_data(DictLoader(thematic))
    aligner.load_reference_data(
        GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
    )

    relevant_distances = np.round(np.arange(0.0, 20.0001, 0.2), 1).tolist()
    process_result = aligner.process(
        relevant_distances=relevant_distances,
        processing_area=processing_area,
    )
    result = process_result.results["theme_1"][1.0]["result"]

    outside_before = thematic_geom.difference(processing_area)
    outside_after = result.difference(processing_area)
    inside_before = thematic_geom.intersection(processing_area)
    inside_after = result.intersection(processing_area)

    print("=== Processing area ADP polygon example ===")
    print("Relevant distances:", f"{relevant_distances[0]}..{relevant_distances[-1]} step 0.2")
    print("Processing area bounds:", processing_area.bounds)
    print("Input area:", round(thematic_geom.area, 3))
    print("Result area:", round(result.area, 3))
    print(
        "Outside scope unchanged:",
        outside_before.equals_exact(outside_after, tolerance=1e-7),
    )
    print(
        "Inside scope changed:",
        not inside_before.equals_exact(inside_after, tolerance=1e-7),
    )
    print("")
    print("=== Per relevant distance summary ===")
    for rd in relevant_distances:
        rd_res = process_result.results["theme_1"][rd]
        rd_geom = rd_res["result"]
        rd_diff = rd_res.get("result_diff")
        rd_diff_area = 0.0 if rd_diff is None or rd_diff.is_empty else float(rd_diff.area)
        rd_remarks = rd_res.get("properties", {}).get(REMARK_FIELD_NAME, [])
        print(
            f"rd={rd:>4.1f} | area={rd_geom.area:>10.3f} | diff_area={rd_diff_area:>10.3f} | remarks={rd_remarks}"
        )

    print("")
    print("=== Evaluate (best prediction) ===")
    evaluated = aligner.evaluate(
        relevant_distances=relevant_distances,
        processing_area=processing_area,
        full_reference_strategy=FullReferenceStrategy.PREFER_FULL_REFERENCE,
    )
    fc = evaluated.get_results_as_geojson(
        result_type=AlignerResultType.EVALUATED_PREDICTIONS,
        add_metadata=True,
        add_original_attributes=True,
        aligner=aligner,
    )
    feat = fc["result"]["features"][0]
    props = feat["properties"]
    print("evaluation:", props.get(EVALUATION_FIELD_NAME))
    print("prediction_score:", props.get(PREDICTION_SCORE))
    print("chosen_rd:", props.get(RELEVANT_DISTANCE_FIELD_NAME))
