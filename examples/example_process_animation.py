import numpy as np

from brdr.aligner import Aligner
from brdr.be.grb.enums import GRBType
from brdr.be.grb.loader import GRBActualLoader
from brdr.loader import GeoJsonLoader
from brdr.viz import animated_map

if __name__ == "__main__":
    # Example: process a range of relevant distances and export animation.
    aligner = Aligner()

    geometry = {  # noqa: F841 - kept for docs: alternative geometry candidate
        "coordinates": [
            [
                [171136.16765265574, 170605.6084498393],
                [171164.53923982676, 170663.78453262427],
                [171071.40019103308, 170714.22290981715],
                [171037.87013346734, 170676.6807086111],
                [171136.16765265574, 170605.6084498393],
            ]
        ],
        "type": "Polygon",
    }

    geometry = {  # Active geometry used in this example.
        "coordinates": [
            [
                [171132.08386359326, 170804.71107799126],
                [171168.04986551206, 170776.33949082025],
                [171144.2635853586, 170743.09601494312],
                [171118.32794254067, 170761.58053385757],
                [171111.52162743648, 170766.9897632298],
                [171125.34919391124, 170785.6892184107],
                [171126.13729355487, 170786.83554516506],
                [171121.9102136481, 170789.8804756064],
                [171132.08386359326, 170804.71107799126],
            ]
        ],
        "type": "Polygon",
    }
    thematic_json = {
        "type": "FeatureCollection",
        "name": "test",
        "crs": {"type": "name", "properties": {"name": "urn:ogc:def:crs:EPSG::31370"}},
        "features": [
            {
                "type": "Feature",
                "properties": {"fid": 1100, "id": 1100, "theme_identifier": "1100"},
                "geometry": geometry,
            }
        ],
    }

    aligner.load_thematic_data(
        GeoJsonLoader(_input=thematic_json, id_property="theme_identifier")
    )
    aligner.load_reference_data(
        GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
    )

    relevant_distances = np.arange(0, 710, 10, dtype=int) / 100
    aligner_result = aligner.process(relevant_distances=relevant_distances)
    thematic_geometries = {
        key: feat.geometry for key, feat in aligner.thematic_data.features.items()
    }
    reference_geometries = {
        key: feat.geometry for key, feat in aligner.reference_data.features.items()
    }

    filename = "animation.gif"
    animated_map(
        aligner_result.results,
        thematic_geometries,
        reference_geometries,
        relevant_distances[-1],
        25,
        150,
        filename,
    )
