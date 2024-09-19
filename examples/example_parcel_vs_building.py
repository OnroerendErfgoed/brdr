import numpy as np

from brdr.aligner import Aligner
from brdr.enums import GRBType
from brdr.grb import GRBActualLoader
from brdr.loader import GeoJsonFileLoader
from brdr.utils import diffs_from_dict_series
from examples import plot_series

# example to check if we can notice if it is better to align to a building instead of a
# parcel
if __name__ == "__main__":
    # Initiate brdr
    aligner_x = Aligner()
    # Load thematic data & reference data (parcels)
    aligner_x.load_thematic_data(
        GeoJsonFileLoader(
            "../tests/testdata/test_parcel_vs_building.geojson", "theme_id"
        )
    )
    aligner_x.load_reference_data(
        GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner_x)
    )  # gebruik de actuele adp-percelen adp= administratieve percelen

    aligner_y = Aligner()
    # Load thematic data & reference data (buildings)
    aligner_y.load_thematic_data(
        GeoJsonFileLoader(
            "../tests/testdata/test_parcel_vs_building.geojson", "theme_id"
        )
    )
    aligner_y.load_reference_data(
        GRBActualLoader(grb_type=GRBType.GBG, partition=1000, aligner=aligner_y)
    )  # gebruik de actuele adp-percelen adp= administratieve percelen

    # Example how to use a series (for histogram)
    series = np.arange(0, 300, 10, dtype=int) / 100
    x_dict_series = aligner_x.process(series, 4, 50)
    x_resulting_areas = diffs_from_dict_series(x_dict_series, aligner_x.dict_thematic)
    y_dict_series = aligner_y.process(series, 4, 50)
    y_resulting_areas = diffs_from_dict_series(y_dict_series, aligner_y.dict_thematic)
    # plot_diffs(series,x_resulting_areas)
    # plot_diffs(series,y_resulting_areas)

    # Make a 1-by-1 comparison for each thematic feature compared to the 2 references (
    # parcels and buildings)
    for key in x_resulting_areas:
        dict_diff = {
            "x" + str(key): x_resulting_areas[key],
            "y" + str(key): y_resulting_areas[key],
        }
        plot_series(series, dict_diff)
