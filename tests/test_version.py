import os

import pytest
import toml

from brdr import __version__

here = os.path.abspath(os.path.dirname(__file__))


@pytest.fixture
def pyproject_version():
    """Fixture to extract the version from pyproject.toml."""
    with open(f"{here}/../pyproject.toml", "r") as f:
        pyproject_data = toml.load(f)
        return pyproject_data["project"]["version"]


def test_version_consistency(pyproject_version):
    """Test to ensure version in pyproject.toml matches __version__."""
    assert pyproject_version == __version__
