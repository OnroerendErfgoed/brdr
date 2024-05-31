import numpy as np
from brdr.aligner import Aligner
from examples import plot_series

# example to check if we can notice if it is better to align to a building instead of a
# parcel
if __name__ == "__main__":
    # Initiate brdr
    x = Aligner()
    # Load thematic data & reference data (parcels)
    x.load_thematic_data_file(
        "../tests/testdata/test_parcel_vs_building.geojson", "theme_id"
    )
    x.load_reference_data_grb_actual(
        grb_type="adp", partition=1000
    )  # gebruik de actuele adp-percelen adp= administratieve percelen

    y = Aligner()
    # Load thematic data & reference data (buildings)
    y.load_thematic_data_file(
        "../tests/testdata/test_parcel_vs_building.geojson", "theme_id"
    )
    y.load_reference_data_grb_actual(
        grb_type="gbg", partition=1000
    )  # gebruik de actuele adp-percelen adp= administratieve percelen

    # Example how to use a series (for histogram)
    series = np.arange(0.1, 5.05, 0.1, dtype=float)
    x_resulting_areas = x.process_series(series, 4, 50)
    y_resulting_areas = y.process_series(series, 4, 50)
    # plot_diffs(series,x_resulting_areas)
    # plot_diffs(series,y_resulting_areas)

    # Make a 1-by-1 comparison for each thematic feature compared to the 2 references (
    # parcels and buildings)
    for key in x_resulting_areas:
        dict_diff = {"x" + str(key): x_resulting_areas[key], "y" + str(key): y_resulting_areas[key]}
        plot_series(series, dict_diff)
