import json
import os

import pytest
import responses

here = os.path.dirname(os.path.abspath(__file__))


@pytest.fixture(autouse=True, scope="session")
def requests_mock():
    """Block all requests calls."""
    with responses.RequestsMock() as rsps:
        rsps.add_passthru("https://www.mercator.vlaanderen.be")
        rsps.add_passthru("https://geo.api.vlaanderen.be")
        rsps.add_passthru("https://inventaris.onroerenderfgoed.be")
        yield rsps


@pytest.fixture(scope="session")
def haspengouw_geojson():
    with open(f"{here}/testdata/haspengouw.geojson", "r") as f:
        return json.loads(f.read())
