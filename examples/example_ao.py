import numpy as np

from brdr.aligner import Aligner
from brdr.enums import GRBType
from brdr.grb import GRBActualLoader
from brdr.oe import OnroerendErfgoedLoader
from examples import show_map, plot_series

if __name__ == "__main__":
    # EXAMPLE to test the algorithm for erfgoedobject with relevant distance 0.2m and
    # od_strategy SNAP_ALL_SIDE

    # Initiate brdr
    aligner = Aligner()
    # Load thematic data & reference data
    aanduidingsobjecten = [117798, 116800, 117881]

    loader = OnroerendErfgoedLoader(aanduidingsobjecten)
    aligner.load_thematic_data(loader)
    loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
    aligner.load_reference_data(loader)

    series = np.arange(0, 500, 20, dtype=int) / 100
    # predict which relevant distances are interesting to propose as resulting geometry
    dict_series, dict_predictions, diffs = aligner.predictor(
        relevant_distances=series, od_strategy=2, threshold_overlap_percentage=50
    )
    for key in dict_predictions.keys():
        diff = {key: diffs[key]}
        plot_series(series, diff)
        show_map(
            {key: dict_predictions[key]},
            {key: aligner.dict_thematic[key]},
            aligner.dict_reference,
        )
