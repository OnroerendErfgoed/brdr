from brdr.aligner import Aligner
from brdr.be.grb.grb import update_featurecollection_to_actual_grb
from brdr.be.grb.loader import GRBFiscalParcelLoader
from brdr.constants import EVALUATION_FIELD_NAME, OBSERVATION_FIELD_NAME
from brdr.loader import GeoJsonLoader

if __name__ == "__main__":
    """
    EXAMPLE of the use of GRB (flanders-specific) function: Update_to_actual_grb
    """
    # Create a featurecollection (aligned on 2022), to use for the 'update_to_actual_grb'
    base_year = "2022"
    base_aligner = Aligner()
    name_thematic_id = "ID"
    loader = GeoJsonLoader(
        _input={
            "type": "FeatureCollection",
            "name": "extract",
            "crs": {
                "type": "name",
                "properties": {"name": "urn:ogc:def:crs:EPSG::31370"},
            },
            "features": [
                {
                    "type": "Feature",
                    "properties": {
                        "testattribute": "test",
                        "nr_calculations": 1,
                        "ID": "206285",
                        "relevant_distance": 2.0,
                        "area": 503.67736346047076,
                        "perimeter": 125.74541473322422,
                        "shape_index": 0.24965468741597097,
                    },
                    "geometry": {
                        "type": "MultiPolygon",
                        "coordinates": [
                            [
                                [
                                    [138539.326299999986077, 193994.138199999986682],
                                    [138529.3663, 193995.566400000010617],
                                    [138522.0997, 193996.6084],
                                    [138514.984399999986636, 193997.6287],
                                    [138505.8261, 193996.615],
                                    [138498.8406, 193996.4314],
                                    [138492.9442, 193996.289500000013504],
                                    [138491.224599999986822, 193996.2481],
                                    [138491.4111, 194004.814699999988079],
                                    [138514.368500000011409, 194005.1297],
                                    [138520.2585, 194004.5753],
                                    [138520.3946, 194005.5833],
                                    [138520.542599999986123, 194009.731999999989057],
                                    [138541.4173, 194007.7292],
                                    [138539.326299999986077, 193994.138199999986682],
                                ]
                            ]
                        ],
                    },
                }
            ],
        },
        id_property="ID",
    )
    base_aligner.load_thematic_data(loader)
    base_aligner.load_reference_data(
        GRBFiscalParcelLoader(year=base_year, aligner=base_aligner)
    )
    aligner_result_base = base_aligner.process(relevant_distances=[2])
    fcs = aligner_result_base.get_results_as_geojson(
        aligner=base_aligner, add_metadata=True, add_original_attributes=True
    )
    featurecollection_base_result = fcs["result"]
    print(featurecollection_base_result)
    # Update Featurecollection to actual version
    featurecollection = update_featurecollection_to_actual_grb(
        featurecollection_base_result,
        id_theme_fieldname="ID",
        base_metadata_field=OBSERVATION_FIELD_NAME,
        max_distance_for_actualisation=3,
    )
    if len(featurecollection) == 0:
        print("empty featurecolection, no updates")
    else:
        # Print results
        for feature in featurecollection["result"]["features"]:
            print(
                feature["properties"]["brdr_id"]
                + ": "
                + feature["properties"][EVALUATION_FIELD_NAME]
            )
        geojson = featurecollection["result"]
        print(geojson)

    featurecollection = update_featurecollection_to_actual_grb(
        featurecollection_base_result,
        id_theme_fieldname="ID",
        base_metadata_field=OBSERVATION_FIELD_NAME,
        max_distance_for_actualisation=0,
    )
    if len(featurecollection) == 0:
        print("empty featurecolection, no updates")
    else:
        # Print results
        for feature in featurecollection["result"]["features"]:
            print(
                feature["properties"][name_thematic_id]
                + ": "
                + feature["properties"][EVALUATION_FIELD_NAME]
            )
        geojson = featurecollection["result"]
        print(geojson)

    featurecollection = update_featurecollection_to_actual_grb(
        featurecollection_base_result,
        id_theme_fieldname="ID",
        max_distance_for_actualisation=3,
    )
    if len(featurecollection) == 0:
        print("empty featurecolection, no updates")
    else:
        # Print results
        for feature in featurecollection["result"]["features"]:
            print(
                feature["properties"][name_thematic_id]
                + ": "
                + feature["properties"][EVALUATION_FIELD_NAME]
            )
        geojson = featurecollection["result"]
        print(geojson)
