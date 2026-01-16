import numpy as np

from brdr.aligner import Aligner
from brdr.be.grb.enums import GRBType
from brdr.be.grb.loader import GRBActualLoader
from brdr.enums import AlignerResultType
from brdr.loader import GeoJsonLoader
from brdr.viz import show_map, plot_difference_by_relevant_distance

# Press the green button in the gutter to run the script.
if __name__ == "__main__":
    """
    EXAMPLE to use the predictor-function to automatically predict which resulting
    geometries are interesting to look at (based on detection of breakpoints and
    relevant distances of 'no-change')
    """

    input_geojson = {
        "type": "FeatureCollection",
        "name": "test_wanted_changes",
        "crs": {"type": "name", "properties": {"name": "urn:ogc:def:crs:EPSG::31370"}},
        "features": [
            {
                "type": "Feature",
                "properties": {"id": 1, "theme_id": "1"},
                "geometry": {
                    "type": "MultiPolygon",
                    "coordinates": [
                        [
                            [
                                [161658.784979313611984, 196026.640970099717379],
                                [161653.250003308057785, 196011.45697008818388],
                                [161615.133042572706472, 196023.607370208803331],
                                [161616.975955285131931, 196029.205066103488207],
                                [161619.898451283574104, 196037.769162107259035],
                                [161620.363155283033848, 196039.130890108644962],
                                [161620.517940590012586, 196039.567879189708037],
                                [161633.552001053554704, 196038.117011958122021],
                                [161631.926488338969648, 196035.414327338279691],
                                [161658.784979313611984, 196026.640970099717379],
                            ]
                        ]
                    ],
                },
            }
        ],
    }
    # Initiate an Aligner
    aligner = Aligner()
    # Load thematic data & reference data
    loader = GeoJsonLoader(_input=input_geojson, id_property="theme_id")

    aligner.load_thematic_data(loader)
    # Load reference data
    loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
    aligner.load_reference_data(loader)

    # PREDICT the 'stable' relevant distances, for a series of relevant distances
    series = np.arange(0, 310, 10, dtype=int) / 100
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
