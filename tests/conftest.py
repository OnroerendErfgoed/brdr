import json
import os
import re
from copy import copy

import pytest
import responses

from tests.testdata.responses import grb_responses, osm_responses
from tests.testdata.responses import inventaris_responses

here = os.path.dirname(os.path.abspath(__file__))


@pytest.fixture(autouse=True, scope="session")
def requests_mock():
    """Block all requests calls."""
    with responses.RequestsMock() as rsps:
        yield rsps


@pytest.fixture(scope="session")
def haspengouw_geojson():
    with open(f"{here}/testdata/haspengouw.geojson", "r") as f:
        return json.loads(f.read())


def multi_url_pattern(*urls):
    pattern = "|".join(f"(^{re.escape(url)})" for url in urls)
    return re.compile(f"{pattern}")


@pytest.fixture
def callback_grb_response(requests_mock):
    response = copy(grb_responses.grb_response)

    def callback(request):
        json_response = json.dumps(response)
        return 200, {}, json_response

    requests_mock.add_callback(
        method=requests_mock.GET,
        url=multi_url_pattern(
            "https://geo.api.vlaanderen.be",
        ),
        callback=callback,
    )
    return response


@pytest.fixture
def mock_grb_response2(callback_grb_response):
    callback_grb_response.update(grb_responses.grb_response2)


@pytest.fixture
def mock_grb_response3(callback_grb_response):
    callback_grb_response.update(grb_responses.grb_response3)


@pytest.fixture
def mock_inventaris_responses(requests_mock):
    requests_mock.add(
        requests_mock.GET,
        multi_url_pattern(
            "https://inventaris.onroerenderfgoed.be/aanduidingsobjecten/131635"
        ),
        json=inventaris_responses.response_131635,
        status=200,
        content_type="application/json",
    )


@pytest.fixture
def mock_be_cadastral_responses(requests_mock):
    with open(f"{here}/testdata/be_cadastral.xml", "r") as f:
        body = f.read()
    requests_mock.add(
        requests_mock.GET,
        multi_url_pattern(
            "https://ccff02.minfin.fgov.be/geoservices/arcgis/services/WMS/Cadastral_LayersWFS/MapServer/WFSServer?"
            # SERVICE=WFS&REQUEST=GetFeature&VERSION=2.0.0&TYPENAMES=CL%3ACadastral_parcel&SRSNAME=http%3A%2F%2Fwww.opengis.net%2Fdef%2Fcrs%2FEPSG%2F0%2F3812&BBOX=672847.835033453%2C671446.8842659479%2C672933.7305729307%2C671548.5887240643%2Chttp%3A%2F%2Fwww.opengis.net%2Fdef%2Fcrs%2FEPSG%2F0%2F3812"
        ),
        body=body,
        status=200,
        content_type="application/xml",
    )


@pytest.fixture
def mock_grb_wfs_capabilities_responses(requests_mock):
    with open(f"{here}/testdata/grb_wfs_capabilities.xml", "r") as f:
        body = f.read()
    requests_mock.add(
        requests_mock.GET,
        multi_url_pattern(
            "https://geo.api.vlaanderen.be/GRB/wfs?service=WFS&version=2.0.0&request=GetCapabilities"
        ),
        body=body,
        status=200,
        content_type="application/xml",
    )


@pytest.fixture
def mock_grb_wfs_responses(requests_mock):
    with open(f"{here}/testdata/grb_wfs_response.xml", "r") as f:
        body = f.read()
    requests_mock.add(
        requests_mock.GET,
        multi_url_pattern(
            "https://geo.api.vlaanderen.be/GRB/wfs?SERVICE=WFS&REQUEST=GetFeature"
        ),
        body=body,
        status=200,
        content_type="application/xml",
    )


@pytest.fixture
def mock_osm_responses(requests_mock):
    requests_mock.add(
        requests_mock.POST,
        "https://overpass-api.de/api/interpreter",
        json=osm_responses.osm_buildings,
        status=200,
        content_type="application/json",
    )
