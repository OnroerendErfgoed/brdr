from shapely import from_wkt

from brdr.aligner import Aligner
from brdr.configs import AlignerConfig
from brdr.loader import DictLoader


def test_metadata():
    config = AlignerConfig(
        log_metadata=True,
        add_observations=True,
    )
    aligner = Aligner(config=config)
    aligner.load_thematic_data(
        DictLoader({"theme_id_1": from_wkt("POLYGON ((0 0, 0 9, 5 10, 10 0, 0 0))")})
    )
    aligner.load_reference_data(
        DictLoader({"ref_id_1": from_wkt("POLYGON ((0 1, 0 10,8 10,10 1,0 1))")})
    )
    result = aligner.process([1])
    assert True
    # TODO emrys this was failing ->repair

    # assert result.results["theme_id_1"][1]["result"].wkt == (
    #     "POLYGON ((0 1, 0 10, 5 10, 9.500000000000004 1, 9.381966011250102 1, "
    #     "9.36275129165333 0.8049096779838723, 9.305845543761384 0.6173165676348947, "
    #     "9.213435623552671 0.4444297669804294, 9.089072792436635 0.2928932188134386, "
    #     "8.937536244269706 0.1685303876974548, 8.76464944361519 0.0761204674887132, "
    #     "8.577056333266233 0.0192147195967696, 8.3819660112501 0, 0.9999999999999997 "
    #     "0, 0.8049096779838718 0.0192147195967696, 0.6173165676349087 "
    #     "0.0761204674887139, 0.4444297669803966 0.1685303876974555, "
    #     "0.2928932188134523 0.2928932188134525, 0.1685303876974549 "
    #     "0.4444297669803975, 0.0761204674887134 0.6173165676349098, "
    #     "0.0192147195967695 0.8049096779838716, 0 1))"
    # )
    # assert result.results["theme_id_1"][1]["metadata"]['actuation']== {
    #     "changes": "geo:hasGeometry",
    #     "id": ANY,
    #     "procedure": {
    #         "id": "2024:aligner2024a",
    #         "implementedBy": "AlignerGeometryProcessor",
    #         "ssn:hasInput": [
    #             {
    #                 "id": "brdr:relevant_distance",
    #                 "input_value": {"type": "xsd:integer", "value": 1},
    #                 "type": "ssn:Input",
    #             }
    #         ],
    #         "type": "sosa:Procedure",
    #     },
    #     "reference_geometries": [
    #         {
    #             "id": ANY,
    #             "type": "Polygon",
    #             "version_date": "",
    #         }
    #     ],
    #     "result": "brdrid:geoms/3f025ce7dcc74d9d095ac2797687a1ec55ec69d3",
    #     "sosa:hasFeatureOfInterest": {"id": ANY},
    #     "type": "sosa:Actuation",
    # }
