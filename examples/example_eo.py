import numpy as np

from brdr.aligner import Aligner
from brdr.enums import GRBType, AlignerResultType
from brdr.grb import GRBActualLoader
from brdr.oe import OnroerendErfgoedLoader, OEType
from brdr.utils import write_geojson
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
    loader = OnroerendErfgoedLoader(objectids=erfgoedobjecten, oetype=OEType.EO)
    aligner.load_thematic_data(loader)
    aligner.load_reference_data(GRBActualLoader(aligner=aligner, grb_type=GRBType.ADP))

    series = np.arange(0, 200, 20, dtype=int) / 100
    # predict which relevant distances are interesting to propose as resulting geometry
    dict_series, dict_predictions, diffs = aligner.predictor(
        relevant_distances=series, od_strategy=2, threshold_overlap_percentage=50
    )
    fcs = aligner.get_results_as_geojson(resulttype=AlignerResultType.PREDICTIONS)
    write_geojson("output/predicted.geojson", fcs["result"])
    write_geojson("output/predicted_diff.geojson", fcs["result_diff"])

    for key in dict_predictions.keys():
        diff = {key: diffs[key]}
        plot_series(series, diff)
        show_map(
            {key: dict_predictions[key]},
            {key: aligner.dict_thematic[key]},
            aligner.dict_reference,
        )
