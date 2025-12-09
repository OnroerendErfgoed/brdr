import pytest

from brdr.aligner import Aligner
from brdr.be.oe.enums import OEType
from brdr.be.oe.loader import OnroerendErfgoedLoader
from tests.conftest import multi_url_pattern
from tests.testdata.responses import mercator_responses


class TestOE:
    def test_onroerenderfgoedloader_by_aanduidid(self, requests_mock):
        requests_mock.add(
            requests_mock.GET,
            "https://www.mercator.vlaanderen.be/raadpleegdienstenmercatorpubliek/wfs?SERVICE=WFS&VERSION=2.0.0&REQUEST=GetFeature&TYPENAMES=ps%3Aps_aandobj&SRSNAME=http%3A%2F%2Fwww.opengis.net%2Fdef%2Fcrs%2FEPSG%2F0%2F31370&outputFormat=application%2Fjson&limit=10000&CQL_FILTER=aanduid_id+IN+%28120288%2C+10275%29",
            json=mercator_responses.response2,
            status=200,
            content_type="application/json",
        )
        loader = OnroerendErfgoedLoader(objectids=[120288, 10275], oetype=OEType.AO)
        aligner = Aligner()
        aligner.load_thematic_data(loader)
        assert len(aligner.dict_thematic.keys()) == 2

    def test_onroerenderfgoedloader_by_erfgoedid(self, requests_mock):
        requests_mock.add(
            requests_mock.GET,
            "https://www.mercator.vlaanderen.be/raadpleegdienstenmercatorpubliek/wfs?SERVICE=WFS&VERSION=2.0.0&REQUEST=GetFeature&TYPENAMES=lu%3Alu_wet_erfgobj_pub&SRSNAME=http%3A%2F%2Fwww.opengis.net%2Fdef%2Fcrs%2FEPSG%2F0%2F31370&outputFormat=application%2Fjson&limit=10000&CQL_FILTER=erfgoed_id+IN+%2842549%29",
            json=mercator_responses.response3,
            status=200,
            content_type="application/json",
        )
        loader = OnroerendErfgoedLoader(objectids=[42549], oetype=OEType.EO)
        aligner = Aligner()
        aligner.load_thematic_data(loader)
        assert len(aligner.dict_thematic.keys()) == 1

    def test_onroerenderfgoedloader_by_bbox(self, requests_mock):
        requests_mock.add(
            requests_mock.GET,
            multi_url_pattern(
                "https://www.mercator.vlaanderen.be/raadpleegdienstenmercatorpubliek/wfs",
            ),
            status=200,
            json=mercator_responses.response4,
            content_type="application/json",
        )
        loader = OnroerendErfgoedLoader(
            bbox=[172000, 172000, 174000, 174000], oetype=OEType.EO
        )
        aligner = Aligner()
        aligner.load_thematic_data(loader)
        assert len(aligner.dict_thematic.keys()) > 0

    def test_onroerenderfgoedloader_by_bbox_and_objectid(self, requests_mock):
        with pytest.raises(Exception):
            OnroerendErfgoedLoader(
                objectids=[42549],
                bbox=[172000, 172000, 174000, 174000],
                oetype=OEType.EO,
            )

        with pytest.raises(Exception):
            OnroerendErfgoedLoader(objectids=None, bbox=None, oetype=OEType.EO)
