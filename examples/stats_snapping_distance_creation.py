import csv
from datetime import datetime
from pathlib import Path

import numpy as np

from brdr.aligner import Aligner
from brdr.loader import GeoJsonFileLoader


if __name__ == "__main__":
    """
    Create a CSV with snapping statistics over a range of relevant distances.
    """
    base_dir = Path(__file__).resolve().parent
    testdata_dir = base_dir.parent / "tests" / "testdata"
    output_dir = base_dir.parent / "tests" / "output"
    output_dir.mkdir(parents=True, exist_ok=True)

    relevant_distances = np.arange(0, 510, 10, dtype=int) / 100

    aligner = Aligner()
    aligner.load_thematic_data(
        GeoJsonFileLoader(str(testdata_dir / "theme_leuven.geojson"), "aanduid_id")
    )
    aligner.load_reference_data(
        GeoJsonFileLoader(str(testdata_dir / "reference_leuven.geojson"), "capakey")
    )

    result = aligner.process(relevant_distances=relevant_distances)
    process_results = result.get_results(aligner=aligner)
    metrics = aligner.get_difference_metrics_for_thematic_data(
        dict_processresults=process_results, thematic_data=aligner.thematic_data
    )

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_path = output_dir / f"stats_{timestamp}.csv"
    with open(output_path, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.writer(csvfile, delimiter=";")
        writer.writerow(
            [
                "distance",
                "key",
                "area",
                "diff_plus",
                "diff_min",
                "diff",
            ]
        )

        for key, by_distance in process_results.items():
            for distance, process_result in by_distance.items():
                if process_result is None:
                    continue
                area = float(process_result["result"].area)
                diff_plus = float(process_result["result_diff_plus"].area)
                diff_min = float(process_result["result_diff_min"].area)
                diff = float(metrics[key][distance])
                writer.writerow([distance, key, area, diff_plus, diff_min, diff])

    print(f"written: {output_path}")
