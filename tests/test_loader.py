from brdr.aligner import Aligner
from brdr.enums import GRBType
from brdr.grb import GRBActualLoader
from brdr.loader import DictLoader, GeoJsonUrlLoader
from brdr.oe import get_oe_dict_by_ids


class TestExamples:

    def test_load_data(self):
        # EXAMPLE
        aligner = Aligner()

        dict_theme = get_oe_dict_by_ids([131635])
        thematic_loader = DictLoader(data_dict=dict_theme)
        reference_loader = GRBActualLoader(
            grb_type=GRBType.ADP, aligner=aligner, partition=1000
        )

        aligner.dict_thematic, props_thematic, thematic_source = (
            thematic_loader.load_data()
        )
        aligner.dict_reference, props_reference, reference_source = (
            reference_loader.load_data()
        )
        assert aligner.dict_reference is not None


def test_load_thematic_data_by_url(requests_mock, haspengouw_geojson):
    requests_mock.add(
        requests_mock.GET,
        "https://mock.com/haspengouw.geojson",
        json=haspengouw_geojson,
        status=200,
    )
    aligner = Aligner()
    aligner.load_thematic_data(
        GeoJsonUrlLoader("https://mock.com/haspengouw.geojson", "Id")
    )
    assert aligner.dict_thematic is not None


def test_load_reference_data_url(requests_mock, haspengouw_geojson):
    requests_mock.add(
        requests_mock.GET,
        "https://mock.com/haspengouw.geojson",
        json=haspengouw_geojson,
        status=200,
    )
    aligner = Aligner()

    aligner.load_reference_data(
        GeoJsonUrlLoader("https://mock.com/haspengouw.geojson", "Id")
    )

    assert aligner.dict_reference is not None
