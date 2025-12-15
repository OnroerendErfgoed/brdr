from unittest.mock import ANY

import numpy as np
import pytest
from shapely import from_wkt

from brdr.aligner import Aligner
from brdr.be.grb.enums import GRBType
from brdr.be.grb.loader import GRBActualLoader
from brdr.enums import FullReferenceStrategy
from brdr.loader import DictLoader


def test_metadata():
    aligner = Aligner()
    aligner.load_thematic_data(
        DictLoader({"theme_id_1": from_wkt("POLYGON ((0 0, 0 9, 5 10, 10 0, 0 0))")})
    )
    aligner.load_reference_data(
        DictLoader({"ref_id_1": from_wkt("POLYGON ((0 1, 0 10,8 10,10 1,0 1))")})
    )
    result = aligner.process([1])

    assert (
        result.results["theme_id_1"][1]["result"].wkt
        == "POLYGON ((0 1, 0 10, 5 10, 9.500000000000004 1, 9.381966011250102 1, 9.36275129165333 0.8049096779838723, 9.305845543761384 0.6173165676348947, 9.213435623552671 0.4444297669804294, 9.089072792436633 0.2928932188134385, 8.937536244269703 0.1685303876974522, 8.764649443615193 0.0761204674887144, 8.577056333266233 0.0192147195967696, 8.3819660112501 0, 0.9999999999999997 0, 0.8049096779838718 0.0192147195967696, 0.6173165676349087 0.0761204674887139, 0.4444297669803966 0.1685303876974555, 0.2928932188134523 0.2928932188134525, 0.1685303876974549 0.4444297669803975, 0.0761204674887134 0.6173165676349098, 0.0192147195967695 0.8049096779838716, 0 1))"
    )
    #TODO test fails because id 'id': 'brdrid:actuations/ is changing at every run
    #assert result.results["theme_id_1"][1]["metadata"] == {'changes': 'geo:hasGeometry', 'id': 'brdrid:actuations/1f72eb49ba2042c481a71e6826777a4e', 'procedure': {'id': '2024:aligner2024a', 'implementedBy': 'AlignerGeometryProcessor', 'ssn:hasInput': [{'id': 'brdr:relevante_afstand', 'input_value': {'type': 'xsd:integer', 'value': 1}, 'type': 'ssn:Input'}], 'type': 'sosa:Procedure'}, 'reference_geometries': [{'id': 'db5753dd9ee041e884265f777cb04395', 'type': 'Polygon', 'version_date': ''}], 'result': 'brdrid:geoms/4143ccf79fed475a00e4294f2d3218dd1e874bcf338cf8fa7e8d802cf439b077', 'sosa:hasFeatureOfInterest': {'id': 'ec9a539b47d145efb115a2da5cf176b3'}, 'type': 'sosa:Actuation'}
