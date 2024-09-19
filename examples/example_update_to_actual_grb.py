from brdr.aligner import Aligner
from brdr.constants import EVALUATION_FIELD_NAME
from brdr.grb import GRBFiscalParcelLoader
from brdr.grb import update_to_actual_grb
from brdr.loader import GeoJsonLoader

# Create a featurecollection (aligned on 2022), to use for the 'update_to_actual_grb'
base_year = "2022"
base_aligner = Aligner()
name_thematic_id = "theme_identifier"
loader = GeoJsonLoader(
    _input={
        "type": "FeatureCollection",
        "name": "extract",
        "crs": {"type": "name", "properties": {"name": "urn:ogc:def:crs:EPSG::31370"}},
        "features": [
            {
                "type": "Feature",
                "properties": {
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
base_process_result = base_aligner.process(relevant_distance=2)
fcs = base_aligner.get_results_as_geojson(formula=True)
featurecollection_base_result = fcs["result"]
print(featurecollection_base_result)
# Update Featurecollection to actual version
featurecollection = update_to_actual_grb(
    featurecollection_base_result, base_aligner.name_thematic_id
)
# Print results
for feature in featurecollection["result"]["features"]:
    print(
        feature["properties"][name_thematic_id]
        + ": "
        + feature["properties"][EVALUATION_FIELD_NAME]
    )
geojson = featurecollection["result"]
print(geojson)
