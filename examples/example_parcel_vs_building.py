import numpy as np

from brdr.aligner import Aligner
from brdr.be.grb.enums import GRBType
from brdr.be.grb.loader import GRBActualLoader
from brdr.loader import GeoJsonFileLoader
from brdr.utils import get_geometry_difference_metrics_from_processresults
from brdr.viz import plot_difference_by_relevant_distance

if __name__ == "__main__":
    """
    # example to check if we can notice if it is better to align to a building instead of a
    # parcel
    """
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
    series = np.arange(0, 310, 10, dtype=int) / 100
    x_aligner_result = aligner_x.process(relevant_distances=series)
    x_results = x_aligner_result.get_results(aligner=aligner_x)
    x_thematic_geometries = {
        key: feat.geometry for key, feat in aligner_x.thematic_data.features.items()
    }
    x_resulting_areas = aligner_x.get_difference_metrics_for_thematic_data(
        dict_processresults=x_results,
        thematic_data=aligner_x.thematic_data,
    )
    y_aligner_result = aligner_y.process(relevant_distances=series)

    y_results = y_aligner_result.get_results(aligner=aligner_x)

    y_resulting_areas = aligner_y.get_difference_metrics_for_thematic_data(
        dict_processresults=y_results, thematic_data=aligner_y.thematic_data
    )

    # Make a 1-by-1 comparison for each thematic feature compared to the 2 references (
    # parcels and buildings)
    for key in x_results.keys():
        diffs_x = get_geometry_difference_metrics_from_processresults(
            dict_processresult=x_results[key],
            geom_thematic=aligner_x.thematic_data.features.get(key).geometry,
            reference_union=aligner_x.reference_data.union,
        )

        diffs_y = get_geometry_difference_metrics_from_processresults(
            dict_processresult=y_results[key],
            geom_thematic=aligner_y.thematic_data.features.get(key).geometry,
            reference_union=aligner_y.reference_data.union,
        )

        dict_diff = {
            "x" + str(key): diffs_x,
            "y" + str(key): diffs_y,
        }
        plot_difference_by_relevant_distance(dict_diff)
