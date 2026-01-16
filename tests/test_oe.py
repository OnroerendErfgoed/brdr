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
            "https://www.mercator.vlaanderen.be/raadpleegdienstenmercatorpubliek/wfs?SERVICE=WFS&VERSION=2.0.0&REQUEST=GetFeature&TYPENAMES=ps%3Aps_aandobj&SRSNAME=http%3A%2F%2Fwww.opengis.net%2Fdef%2Fcrs%2FEPSG%2F0%2F31370&outputFormat=application%2Fjson&limit=10000&CQL_FILTER=uri+IN+%28%27https%3A%2F%2Fid.erfgoed.net%2Faanduidingsobjecten%2F120288%27%2C%27https%3A%2F%2Fid.erfgoed.net%2Faanduidingsobjecten%2F10275%27%29",
            json=mercator_responses.response2,
            status=200,
            content_type="application/json",
        )
        loader = OnroerendErfgoedLoader(
            objectids=[
                "https://id.erfgoed.net/aanduidingsobjecten/120288",
                "https://id.erfgoed.net/aanduidingsobjecten/10275",
            ],
            oetype=OEType.AO,
        )
        aligner = Aligner()
        aligner.load_thematic_data(loader)
        assert len(aligner.thematic_data.features.keys()) == 2

    def test_onroerenderfgoedloader_by_erfgoedid(self, requests_mock):
        requests_mock.add(
            requests_mock.GET,
            "https://www.mercator.vlaanderen.be/raadpleegdienstenmercatorpubliek/wfs?SERVICE=WFS&VERSION=2.0.0&REQUEST=GetFeature&TYPENAMES=lu:lu_wet_erfgobj_pub&SRSNAME=http://www.opengis.net/def/crs/EPSG/0/31370&outputFormat=application/json&limit=10000&CQL_FILTER=uri+IN+('https://id.erfgoed.net/erfgoedobjecten/42549')",
            json=mercator_responses.response3,
            status=200,
            content_type="application/json",
        )
        loader = OnroerendErfgoedLoader(
            objectids=[
                "https://id.erfgoed.net/erfgoedobjecten/42549",
            ],
            oetype=OEType.EO,
        )
        aligner = Aligner()
        aligner.load_thematic_data(loader)
        assert len(aligner.thematic_data.features.keys()) == 1

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
        assert len(aligner.thematic_data.features.keys()) > 0

    def test_onroerenderfgoedloader_by_bbox_and_objectid(self, requests_mock):
        with pytest.raises(Exception):
            OnroerendErfgoedLoader(
                objectids=[42549],
                bbox=[172000, 172000, 174000, 174000],
                oetype=OEType.EO,
            )

        with pytest.raises(Exception):
            OnroerendErfgoedLoader(objectids=None, bbox=None, oetype=OEType.EO)
