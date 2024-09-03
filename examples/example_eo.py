import numpy as np

from brdr.aligner import Aligner
from brdr.enums import GRBType
from brdr.utils import get_oe_dict_by_ids, write_geojson, dict_series_by_keys
from examples import show_map, plot_series

if __name__ == "__main__":
    # EXAMPLE to test the algorithm for erfgoedobject with relevant distance 0.2m and
    # od_strategy SNAP_ALL_SIDE

    # Initiate brdr
    aligner = Aligner()
    # Load thematic data & reference data
    # dict_theme = get_oe_dict_by_ids([206363], oetype='erfgoedobjecten')

    erfgoedobjecten = [
        # 206407,
        # 206403,
        # 206372,
        # 206369,
        # 206377,
        # 206371,
        # 206370,
        # 206368,
        206786
    ]
    dict_theme = get_oe_dict_by_ids(erfgoedobjecten, oetype="erfgoedobjecten")
    aligner.load_thematic_data_dict(dict_theme)
    aligner.load_reference_data_grb_actual(grb_type=GRBType.ADP, partition=1000)

    # RESULTS
    # rel_dist = 0.2
    # dict_results_by_distance = {}
    # #put resulting tuple in a dictionary
    # dict_results_by_distance[rel_dist] = aligner.process_dict_thematic(rel_dist,2)
    # aligner.export_results("output/")
    # show_map(dict_results_by_distance, aligner.dict_thematic, aligner.dict_reference)

    series = np.arange(0, 200, 20, dtype=int) / 100
    # predict which relevant distances are interesting to propose as resulting geometry
    dict_series, dict_predicted, diffs = aligner.predictor(
        relevant_distances=series, od_strategy=2, threshold_overlap_percentage=50
    )
    fcs = aligner.get_predictions_as_geojson(series_dict=dict_predicted)
    write_geojson("output/predicted.geojson", fcs["result"])
    write_geojson("output/predicted_diff.geojson", fcs["result_diff"])

    dict_predicted = dict_series_by_keys(dict_predicted)
    for key in dict_predicted.keys():
        diff = {key: diffs[key]}
        plot_series(series, diff)
        show_map(
            dict_predicted[key],
            {key: aligner.dict_thematic[key]},
            aligner.dict_reference,
        )
