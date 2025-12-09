import pytest

from brdr.aligner import Aligner
from brdr.be.grb.enums import GRBType
from brdr.be.grb.loader import GRBActualLoader
from brdr.be.oe.loader import OnroerendErfgoedLoader

from brdr.loader import GeoJsonUrlLoader
from tests.testdata.responses import mercator_responses


class TestExamples:

    @pytest.mark.usefixtures("callback_grb_response")
    def test_load_data(self,requests_mock):
        requests_mock.add(
            requests_mock.GET,
            "https://www.mercator.vlaanderen.be/raadpleegdienstenmercatorpubliek/wfs?SERVICE=WFS&VERSION=2.0.0&REQUEST=GetFeature&TYPENAMES=ps%3Aps_aandobj&SRSNAME=http%3A%2F%2Fwww.opengis.net%2Fdef%2Fcrs%2FEPSG%2F0%2F31370&outputFormat=application%2Fjson&limit=10000&CQL_FILTER=aanduid_id+IN+%28131635%29",
            json=mercator_responses.response1,
            status=200,
        )
        # EXAMPLE
        aligner = Aligner()
        thematic_loader = OnroerendErfgoedLoader(objectids=[131635])
        reference_loader = GRBActualLoader(
            grb_type=GRBType.ADP, aligner=aligner, partition=1000
        )

        aligner.dict_thematic, props_thematic, thematic_source = (
            thematic_loader.load_data()
        )
        aligner.dict_reference, props_reference, reference_source = (
            reference_loader.load_data()
        )
        assert aligner.dict_thematic is not None
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
