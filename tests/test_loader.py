import pytest

from brdr.aligner import Aligner
from brdr.be.grb.enums import GRBType
from brdr.be.grb.loader import GRBActualLoader
from brdr.be.oe.utils import get_oe_dict_by_ids
from brdr.loader import DictLoader
from brdr.loader import GeoJsonUrlLoader


class TestExamples:

    @pytest.mark.usefixtures("callback_grb_response")
    @pytest.mark.usefixtures("mock_inventaris_responses")
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
