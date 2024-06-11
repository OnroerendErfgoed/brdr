import numpy as np
from brdr.aligner import Aligner
from brdr.utils import filter_resulting_series_by_key
from examples import show_map

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    """
    example to use the predictor-function to automatically predict which resulting geometries are interesting to look at 
    (based on detection of breakpoints and relevant distances of 'no-change')
    """
    # TODO: future possibilities to choose the best proposal automatically? (AI-machine learning?)
    ##Initiate a Aligner
    aligner = Aligner()
    ##Load thematic data & reference data
    aligner.load_thematic_data_file("../tests/testdata/test_wanted_changes.geojson", 'theme_id')
    aligner.load_reference_data_grb_actual(grb_type='adp', partition=1000) #gebruik de actuele adp-percelen adp= administratieve percelen

    series = np.arange(0.1, 5.00, 0.2, dtype=float)
    #predict which relevant distances are interesting to propose as resulting geometry
    dict_predicted = aligner.predictor(relevant_distances=series, od_strategy=4,treshold_overlap_percentage=50)
    for key in dict_predicted.keys():
        show_map(dict_predicted[key], {key:aligner.dict_thematic[key]}, aligner.dict_reference)




