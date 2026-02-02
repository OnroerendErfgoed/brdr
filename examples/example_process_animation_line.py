from datetime import datetime

import numpy as np

from brdr.aligner import Aligner
from brdr.be.grb.enums import GRBType
from brdr.be.grb.loader import GRBActualLoader
from brdr.loader import GeoJsonLoader
from brdr.viz import animated_map

if __name__ == "__main__":
    # EXAMPLE to process a series of relevant distances
    # Initiate brdr
    aligner = Aligner()
    # Load thematic data
    geometry = {"coordinates":[[[171647.88034171,226384.40524778998],[171649.5703417,226381.8352478],[171650.0003417,226381.1752478],[171726.18034169,226420.87524778],[171748.02034169,226433.39524778],[171769.37034168,226441.33524778],[171826.76034167,226461.27524777]]],"type":"MultiLineString"}

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

    loader = GeoJsonLoader(_input=thematic_json, id_property="theme_identifier")
    aligner.load_thematic_data(loader)
    # Load reference data: The actual GRB-parcels
    aligner.load_reference_data(
        GRBActualLoader(grb_type=GRBType.Wegsegment, partition=1000, aligner=aligner)
    )
    # PROCESS a series of relevant distances
    relevant_distances = np.arange(0, 710, 10, dtype=int) / 100
    aligner_result = aligner.process(
        relevant_distances=relevant_distances,
    )
    thematic_geometries = {
        key: feat.geometry for key, feat in aligner.thematic_data.features.items()
    }
    reference_geometries = {
        key: feat.geometry  for key, feat in aligner.reference_data.features.items()
    }
    # SHOW results: map and plotted changes
    now = datetime.now()  # current date and time
    date_time = now.strftime("%m_%d_%Y_%H_%M_%S")
    # filename = "animation_" + str(date_time) + ".gif"
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
