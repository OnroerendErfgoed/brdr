import numpy as np

from brdr.aligner import Aligner
from examples import plot_diffs
from examples import show_results

# example to test if we can separate wanted deviations from unwanted deviations.
# By calculating results with a series of relevant distances we get a graphic.
# There is a difference in the graphic for wanted and unwanted deviations:
#   *unwanted deviations: will evaluate more gradually
#   *wanted deviations: have a harder breakpoint
# The type of breakpoint, breakpoint-distance and length between breakpoints can give
# an indication of wanted and unwanted deviations

# TODO: future possibilities to detect the 'wanted' breakpoint automatically? (
#  AI-machine learning?)


if __name__ == "__main__":
    # Initiate brdr
    x = Aligner()
    # Load thematic data & reference data
    x.load_thematic_data_file(
        "../tests/testdata/test_wanted_changes.geojson", "theme_id"
    )
    x.load_reference_data_grb_actual(
        grb_type="adp", partition=1000
    )  # gebruik de actuele adp-percelen adp= administratieve percelen

    # Example how to use the Aligner
    r, rd, rd_plus, rd_min, sd, si = x.process_dict_thematic(2, 4)
    x.export_results("output/")
    show_results(r, rd_plus,rd_min)

    # Possibility to get the descriptive formula of a thematic feature
    # for key in r:
    #     x.get_formula(r[key])

    # Example how to use a series (for histogram)
    # series = [0.2, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 6, 8, 10]
    # series = [0.1,0.2,0.3,0.4, 0.5,1,2,3,4]
    series = np.arange(0.1, 10.05, 0.1, dtype=float)
    resulting_areas = x.process_series(series, 4, 50)
    plot_diffs(series, resulting_areas)
