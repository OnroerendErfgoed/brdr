import numpy as np
from brdr.aligner import Aligner
from brdr.utils import get_breakpoints_zerostreak, diffs_from_dict_series
from examples import show_results, plot_series, show_individual_results

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    """
    example to test if we can seperate wanted deviations from unwanted deviations.
    By calculating results with a series of relevant distances we get a graphic.
    There is a difference in the graphic for wanted and unwanted deviations:
    * unwanted deviations: will evoluate more gradually
    * wanted deviations: have a harder breakpoint
    The type of breakpoint, breakpoint-distance and length between breakpoints can gave an indication of wanted and unwanted deviations
    """
    # TODO: future possibilities to detect the 'wanted' breakpoint automatically? (AI-machine learning?)
    ##Initiate a Aligner
    aligner = Aligner()
    ##Load thematic data & reference data
    aligner.load_thematic_data_file("../tests/testdata/test_wanted_changes.geojson", 'theme_id')
    aligner.load_reference_data_grb_actual(grb_type='adp', partition=1000) #gebruik de actuele adp-percelen adp= administratieve percelen

    #Example how to use a series (for histogram)
    #series = [0.2, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 6, 8, 10]
    #series = [0.1,0.2,0.3,0.4, 0.5,1,2,3,4]
    series = np.arange(0.1, 5.00, 0.1, dtype=float)
    dict_predicted = aligner.predictor(relevant_distances=series, od_strategy=4,full_overlap_percentage=50)
    for key in dict_predicted.keys():
        for dist in dict_predicted[key]:
            show_individual_results(dict_predicted[key][dist][0],dict_predicted[key][dist][2],dict_predicted[key][dist][3],aligner.dict_thematic[key],aligner.dict_reference)




