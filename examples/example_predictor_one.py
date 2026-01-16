from brdr.aligner import Aligner
from brdr.be.grb.enums import GRBType
from brdr.be.grb.loader import GRBActualLoader
from brdr.enums import AlignerResultType
from brdr.loader import GeoJsonFileLoader
from brdr.viz import show_map, plot_difference_by_relevant_distance

# Press the green button in the gutter to run the script.
if __name__ == "__main__":
    """
    EXAMPLE to use the predictor-function to automatically predict which resulting
    geometries are interesting to look at (based on detection of breakpoints and
    relevant distances of 'no-change')
    """
    # Initiate an Aligner
    aligner = Aligner()
    # Load thematic data & reference data
    loader = GeoJsonFileLoader("input/predictor_one.geojson", "theme_id")

    aligner.load_thematic_data(loader)
    # Load reference data
    loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
    aligner.load_reference_data(loader)

    # PREDICT the 'stable' relevant distances, for a series of relevant distances
    series = [3]
    # predict which relevant distances are interesting to propose as resulting geometry
    aligner_result = aligner.predict(
        relevant_distances=series,
    )
    dict_predictions = aligner_result.get_results(
        aligner=aligner, result_type=AlignerResultType.PREDICTIONS
    )

    # SHOW results of the predictions
    fcs = aligner_result.get_results_as_geojson(add_metadata=False, aligner=aligner)
    diffs_dict = aligner.get_difference_metrics_for_thematic_data(
        dict_processresults=aligner_result.results, thematic_data=aligner.thematic_data
    )
    reference_geometries = {
        key: feat.geometry for key, feat in aligner.reference_data.features.items()
    }
    if fcs is None or "result" not in fcs:
        print("empty predictions")
    else:
        print(fcs["result"])
        for key in dict_predictions:
            plot_difference_by_relevant_distance({key: diffs_dict[key]})
            show_map(
                {key: dict_predictions[key]},
                {key: aligner.thematic_data.features[key].geometry},
                reference_geometries,
            )
