from shapely import affinity
from shapely.geometry import box

from brdr.aligner import Aligner
from brdr.loader import DictLoader


if __name__ == "__main__":
    """
    Example: scope alignment to a selected area with `processing_area`.

    Use case:
    - A thematic geometry spans multiple zones.
    - Only one zone should be aligned now.
    - The remainder must stay exactly as-is.
    """

    # 1) Create synthetic thematic + reference data (no network dependency).
    # Thematic geometry is intentionally shifted to the right by 1 meter.
    thematic_geom_original = affinity.translate(box(0, 0, 20, 10), xoff=1.0, yoff=0.0)
    thematic = {"theme_1": thematic_geom_original}

    # Reference geometry is split in two parts to mimic separate boundaries.
    reference = {
        "ref_left": box(0, 0, 10, 10),
        "ref_right": box(10, 0, 20, 10),
    }

    aligner = Aligner(crs="EPSG:31370")
    aligner.load_thematic_data(DictLoader(data_dict=thematic))
    aligner.load_reference_data(DictLoader(data_dict=reference))

    # 2) Define a scope: only left half should be actively aligned.
    processing_area = box(0, -5, 10, 15)

    rd = 1.5
    full_result = aligner.process(relevant_distances=[rd]).results["theme_1"][rd][
        "result"
    ]
    scoped_result = aligner.process(
        relevant_distances=[rd],
        processing_area=processing_area,
    ).results["theme_1"][rd]["result"]

    # 3) Check left vs right side behavior.
    left_half = box(0, -5, 10, 15)
    right_half = box(10, -5, 25, 15)

    scoped_left = scoped_result.intersection(left_half)
    scoped_right = scoped_result.intersection(right_half)
    original_right = thematic_geom_original.intersection(right_half)

    print("=== Processing area example ===")
    print(f"Relevant distance: {rd}")
    print(f"Processing area bounds: {processing_area.bounds}")
    print("")
    print("Full processing (no scope):")
    print(f"  result bounds: {full_result.bounds}")
    print("")
    print("Scoped processing:")
    print(f"  result bounds: {scoped_result.bounds}")
    print(
        "  right side unchanged vs original:",
        scoped_right.equals_exact(original_right, tolerance=1e-9),
    )
    print(
        "  left side changed vs original:",
        not scoped_left.equals_exact(
            thematic_geom_original.intersection(left_half), tolerance=1e-9
        ),
    )
