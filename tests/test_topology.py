import unittest

from shapely.geometry import Polygon

from brdr.aligner import Aligner
from brdr.enums import GRBType
from brdr.grb import (
    GRBActualLoader,
)
from brdr.loader import GeoJsonLoader


class TestTopology(unittest.TestCase):
    def setUp(self):
        # Create a sample geometry for testing
        self.sample_aligner = Aligner()
        self.sample_geom = Polygon([(0, 0), (0, 10), (10, 10), (10, 0)])

    def test_topology(self):
        """
        Test if parameter preserve_topology is working"""
        # Initiate an Aligner
        geojson = {
            "type": "FeatureCollection",
            "name": "two",
            "crs": {
                "type": "name",
                "properties": {"name": "urn:ogc:def:crs:EPSG::31370"},
            },
            "features": [
                {
                    "type": "Feature",
                    "properties": {
                        "OIDN": 5763414,
                        "UIDN": 15089212,
                        "VERSIE": 1,
                        "BEGINDATUM": "2021-05-16",
                        "VERSDATUM": "2021-05-16",
                        "EINDDATUM": "2024-02-05",
                        "CAPAKEY": "24126B0049/00Z000",
                        "CANU": "49Z",
                        "FISCDATUM": "2022-01-01",
                        "BHRDR": 2,
                        "LBLBHRDR": "AAPD",
                        "NISCODE": "24062",
                        "BGNINV": 11,
                        "LBLBGNINV": "adpupdate",
                        "EINDINV": 11,
                        "LBLEINDINV": "adpupdate",
                        "BEWERK": 1,
                        "LBLBEWERK": "verwijderd",
                        "LENGTE": 117.62,
                        "OPPERVL": 398.32,
                    },
                    "geometry": {
                        "type": "MultiPolygon",
                        "coordinates": [
                            [
                                [
                                    [174193.075163975358009, 179501.127230677753687],
                                    [174193.069019973278046, 179501.100158676505089],
                                    [174183.54901996999979, 179504.310078680515289],
                                    [174184.360603965818882, 179506.715262681245804],
                                    [174192.349019974470139, 179530.390142697840929],
                                    [174197.229019977152348, 179552.460094712674618],
                                    [174199.52905198186636, 179551.690110713243484],
                                    [174204.018971979618073, 179550.170110709965229],
                                    [174196.94549997895956, 179518.199998687952757],
                                    [174196.899035975337029, 179517.990142688155174],
                                    [174193.623707972466946, 179503.546174678951502],
                                    [174193.075163975358009, 179501.127230677753687],
                                ]
                            ]
                        ],
                    },
                },
                {
                    "type": "Feature",
                    "properties": {
                        "OIDN": 3322239,
                        "UIDN": 15109808,
                        "VERSIE": 3,
                        "BEGINDATUM": "2011-12-08",
                        "VERSDATUM": "2021-05-16",
                        "EINDDATUM": "2024-02-06",
                        "CAPAKEY": "24126B0056/00F000",
                        "CANU": "56F",
                        "FISCDATUM": "2022-01-01",
                        "BHRDR": 2,
                        "LBLBHRDR": "AAPD",
                        "NISCODE": "24062",
                        "BGNINV": 11,
                        "LBLBGNINV": "adpupdate",
                        "EINDINV": 11,
                        "LBLEINDINV": "adpupdate",
                        "BEWERK": 3,
                        "LBLBEWERK": "geometriewijziging, niet beduidend",
                        "LENGTE": 290.35,
                        "OPPERVL": 2313.54,
                    },
                    "geometry": {
                        "type": "MultiPolygon",
                        "coordinates": [
                            [
                                [
                                    [174245.362332008779049, 179581.87654273211956],
                                    [174239.231516011059284, 179554.011070713400841],
                                    [174227.245403997600079, 179557.761854715645313],
                                    [174219.442075990140438, 179560.203838717192411],
                                    [174210.462427988648415, 179520.346942689269781],
                                    [174206.717275984585285, 179509.097086682915688],
                                    [174204.78287598490715, 179503.286270678043365],
                                    [174202.223963983356953, 179497.917822673916817],
                                    [174202.209243983030319, 179497.886910676956177],
                                    [174193.069019973278046, 179501.100158676505089],
                                    [174193.075163975358009, 179501.127230677753687],
                                    [174193.623707972466946, 179503.546174678951502],
                                    [174196.899035975337029, 179517.990142688155174],
                                    [174196.94549997895956, 179518.199998687952757],
                                    [174204.018971979618073, 179550.170110709965229],
                                    [174207.168987981975079, 179564.370046719908714],
                                    [174200.35887598246336, 179566.609854724258184],
                                    [174188.40898796916008, 179570.540030725300312],
                                    [174193.813019976019859, 179598.002366743981838],
                                    [174245.362332008779049, 179581.87654273211956],
                                ]
                            ]
                        ],
                    },
                },
            ],
        }

        aligner = Aligner(crs="EPSG:31370", preserve_topology=True)
        loader = GeoJsonLoader(_input=geojson, id_property="CAPAKEY")
        aligner.load_thematic_data(loader)
        aligner.load_reference_data(
            GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
        )
        relevant_distance = 2
        process_result = aligner.process(
            relevant_distance=relevant_distance,
        )

        self.assertEqual(len(process_result), 2)
        dict_predictions_evaluated, prop_dictionary = aligner.evaluate()
        print(dict_predictions_evaluated)
        print(prop_dictionary)
