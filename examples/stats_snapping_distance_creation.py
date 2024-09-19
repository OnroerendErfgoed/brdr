import csv
import logging
import time

import numpy as np

from brdr.aligner import Aligner
from brdr.loader import GeoJsonFileLoader

# Code to create stats.csv

time = str(time.time())
array_threshold_overlap_percentage = [50]
array_od = [1]
# array_relevant_distance = [0.2, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 8, 10]
array_relevant_distance = np.arange(0, 500, 10, dtype=int) / 100
print(array_relevant_distance)

x = Aligner()
x.load_thematic_data(
    GeoJsonFileLoader("../tests/testdata/theme_leuven.geojson", "aanduid_id")
)
x.load_reference_data(
    GeoJsonFileLoader("../tests/testdata/reference_leuven.geojson", "capakey")
)
with open("../tests/output/stats" + time + ".csv", "w", newline="") as csvfile:
    writer = csv.writer(
        csvfile, delimiter=";"
    )  # ,quotechar='"', quoting=csv.QUOTE_MINIMAL)
    writer.writerow(
        [
            "distance",
            "strategy",
            "full_percentage",
            "key",
            "area",
            "diff_plus",
            "diff_min",
            "diff",
            "diff_prc",
        ]
    )

    for full_percentage in array_threshold_overlap_percentage:
        logging.info("full overlap percentage: " + str(full_percentage))
        for od in array_od:
            logging.info("od_strategy: " + str(od))
            for s in array_relevant_distance:
                logging.info("relevant_distance: " + str(s))
                (
                    results,
                    results_diff,
                    results_diff_plus,
                    results_diff_min,
                    si,
                    sd,
                ) = x.process(s, od, full_percentage)
                for key in results:
                    results_area = results[key].area
                    results_diff_plus_area = results_diff_plus[key].area
                    results_diff_min_area = results_diff_min[key].area
                    diff_area = results_diff_plus_area + results_diff_min_area
                    diff_area_prc = diff_area * 100 / results_area
                    writer.writerow(
                        [
                            s,
                            od,
                            full_percentage,
                            key,
                            results_area,
                            results_diff_plus_area,
                            results_diff_min_area,
                            diff_area,
                            diff_area_prc,
                        ]
                    )
