import numpy as np

from brdr.aligner import Aligner
from brdr.be.grb.enums import GRBType
from brdr.be.grb.loader import GRBActualLoader
from brdr.enums import AlignerResultType
from brdr.loader import GeoJsonLoader
from brdr.viz import plot_difference_by_relevant_distance, show_map

if __name__ == "__main__":
    """Example: generate predicted candidate geometries from distance series."""

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

    aligner = Aligner()
    aligner.load_thematic_data(
        GeoJsonLoader(_input=input_geojson, id_property="theme_id")
    )
    aligner.load_reference_data(
        GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
    )

    relevant_distances = np.arange(0, 310, 10, dtype=int) / 100
    aligner_result = aligner.predict(relevant_distances=relevant_distances)
    predictions = aligner_result.get_results(
        aligner=aligner, result_type=AlignerResultType.PREDICTIONS
    )

    geojson_results = aligner_result.get_results_as_geojson(
        add_metadata=False, aligner=aligner
    )
    diffs_dict = aligner.get_difference_metrics_for_thematic_data(
        dict_processresults=aligner_result.results, thematic_data=aligner.thematic_data
    )
    reference_geometries = {
        key: feat.geometry for key, feat in aligner.reference_data.features.items()
    }

    if geojson_results is None or "result" not in geojson_results:
        print("empty predictions")
    else:
        print(geojson_results["result"])
        for key in predictions:
            plot_difference_by_relevant_distance({key: diffs_dict[key]})
            show_map(
                {key: predictions[key]},
                {key: aligner.thematic_data.features[key].geometry},
                reference_geometries,
            )
