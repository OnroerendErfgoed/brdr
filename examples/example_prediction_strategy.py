import numpy as np

from brdr.aligner import Aligner
from brdr.constants import EVALUATION_FIELD_NAME
from brdr.enums import GRBType, FullStrategy, Evaluation
from brdr.grb import GRBActualLoader
from brdr.loader import GeoJsonLoader

thematic_json = {
    "type": "FeatureCollection",
    "name": "uc5",
    "crs": {"type": "name", "properties": {"name": "urn:ogc:def:crs:EPSG::31370"}},
    "features": [
        {
            "type": "Feature",
            "properties": {"fid": 1, "dossiernummer": "D_68365"},
            "geometry": {
                "type": "Polygon",
                "coordinates": [
                    [
                        [172648.411900000093738, 171970.078100000013364],
                        [172641.128199999919161, 171979.2324],
                        [172638.179699999920558, 171982.9382],
                        [172633.854899999918416, 171988.3738],
                        [172633.317000841663685, 171989.049898942059372],
                        [172623.4652, 172001.432],
                        [172610.793499999970663, 172017.358100000128616],
                        [172610.942697853111895, 172017.477398283459479],
                        [172650.691, 172049.2673],
                        [172654.167900000087684, 172052.047999999922467],
                        [172656.0643, 172049.772999999986496],
                        [172657.2538, 172048.346],
                        [172677.110299999912968, 172024.525499999901513],
                        [172677.8891, 172023.591299999912735],
                        [172692.862399999867193, 172005.628700000088429],
                        [172693.374199999903794, 172006.038],
                        [172696.634700000082375, 172001.575399999826914],
                        [172709.757699999987381, 171986.0387],
                        [172698.4481, 171976.4558],
                        [172665.904004495765548, 171948.880303809302859],
                        [172665.794206234189915, 171948.787305280333385],
                        [172665.529999999969732, 171948.563400000159163],
                        [172648.411900000093738, 171970.078100000013364],
                    ]
                ],
            },
        },
        {
            "type": "Feature",
            "properties": {"fid": 2, "dossiernummer": "D_352"},
            "geometry": {
                "type": "Polygon",
                "coordinates": [
                    [
                        [172443.633200000098441, 171905.855499999917811],
                        [172450.522299999924144, 171922.915599999832921],
                        [172471.899299999902723, 171914.2226],
                        [172458.880800000013551, 171881.983800000074552],
                        [172447.773599999636644, 171885.025300000183051],
                        [172436.992100000090431, 171889.409699999901932],
                        [172443.633200000098441, 171905.855499999917811],
                    ]
                ],
            },
        },
        {
            "type": "Feature",
            "properties": {"fid": 3, "dossiernummer": "D_28311"},
            "geometry": {
                "type": "Polygon",
                "coordinates": [
                    [
                        [172179.29800000009709, 171788.713699999905657],
                        [172182.4732, 171780.029],
                        [172157.412400000175694, 171770.505299999989802],
                        [172129.436300000088522, 171848.2427],
                        [172129.052500541059999, 171849.309098496683873],
                        [172128.886999999987893, 171849.769],
                        [172148.2379, 171862.7293],
                        [172151.4504, 171864.880899999901885],
                        [172151.894199230970116, 171863.666902103519533],
                        [172179.29800000009709, 171788.713699999905657],
                    ]
                ],
            },
        },
    ],
}
# 4001 will give a result, 4002 &4009 should also give the original geometry


aligner = Aligner(relevant_distances=np.arange(0, 310, 10, dtype=int) / 100)

loader = GeoJsonLoader(_input=thematic_json, id_property="fid")
aligner.load_thematic_data(loader)
# Load reference data: The actual GRB-parcels
loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
aligner.load_reference_data(loader)

# Use the EVALUATE-function
dict_predictions_evaluated, prop_dictionary = aligner.evaluate(
    max_predictions=1,
    full_strategy=FullStrategy.ONLY_FULL,
    multi_to_best_prediction=True,
)

assert prop_dictionary[1][0][EVALUATION_FIELD_NAME] == Evaluation.PREDICTION_UNIQUE_FULL
assert prop_dictionary[2][0][EVALUATION_FIELD_NAME] == Evaluation.TO_CHECK_NO_PREDICTION
assert prop_dictionary[3][0][EVALUATION_FIELD_NAME] == Evaluation.TO_CHECK_NO_PREDICTION
