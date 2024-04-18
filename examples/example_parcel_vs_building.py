import numpy as np
from matplotlib import pyplot as plt

from brdr.auto_referencer import AutoReferencer

# example to check if we can notice if it is better to align to a building instead of a
# parcel
if __name__ == "__main__":
    # Initiate brdr
    x = AutoReferencer()
    # Load thematic data & reference data (parcels)
    x.load_thematic_data_file(
        "../tests/testdata/test_parcel_vs_building.geojson", "theme_id"
    )
    x.load_reference_data_grb_actual(
        grb_type="adp", partition=1000
    )  # gebruik de actuele adp-percelen adp= administratieve percelen

    y = AutoReferencer()
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
        if len(x_resulting_areas[key]) == len(series):
            lst_diffs = list(x_resulting_areas[key].values())
            plt.plot(series, lst_diffs, label="x" + str(key))
        if len(y_resulting_areas[key]) == len(series):
            lst_diffs = list(y_resulting_areas[key].values())
            plt.plot(series, lst_diffs, label="y" + str(key))
        plt.xlabel("relevante afstand")
        plt.ylabel("diff")
        plt.title("relevante afstand vs diff")
        plt.legend()
        plt.show()
