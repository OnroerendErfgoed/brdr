import numpy as np

from brdr.aligner import Aligner
from brdr.utils import get_oe_dict_by_ids
from examples import show_map, plot_series

if __name__ == "__main__":
    # EXAMPLE to test the algorithm for erfgoedobject with relevant distance 0.2m and od_strategy SNAP_ALL_SIDE

    # Initiate brdr
    aligner = Aligner()
    # Load thematic data & reference data
    #dict_theme = get_oe_dict_by_ids([206363], oetype='erfgoedobjecten')
    aanduidingsobjecten = range(1,10)
    dict_theme = get_oe_dict_by_ids(aanduidingsobjecten, oetype='aanduidingsobjecten')
    aligner.load_thematic_data_dict(dict_theme)
    aligner.load_reference_data_grb_actual(grb_type="adp", partition=1000)

    #RESULTS
    # rel_dist = 0.2
    # dict_results_by_distance = {}
    # #put resulting tuple in a dictionary
    # dict_results_by_distance[rel_dist] = aligner.process_dict_thematic(rel_dist,2)
    # aligner.export_results("output/")
    # show_map(dict_results_by_distance, aligner.dict_thematic, aligner.dict_reference)

    series = np.arange(0, 500, 20, dtype=int)/100
    #predict which relevant distances are interesting to propose as resulting geometry
    dict_predicted, diffs = aligner.predictor(relevant_distances=series, od_strategy=2,treshold_overlap_percentage=50)
    for key in dict_predicted.keys():
        diff ={}
        diff[key]= diffs[key]
        plot_series(series,diff)
        show_map(dict_predicted[key], {key:aligner.dict_thematic[key]}, aligner.dict_reference)