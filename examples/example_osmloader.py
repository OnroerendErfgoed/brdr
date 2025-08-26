from brdr.aligner import Aligner
from brdr.geometry_utils import geom_from_wkt
from brdr.loader import DictLoader
from brdr.osm import OSMLoader

if __name__ == "__main__":
    """
    Example to reference data with OSMLoader
    """
    # CREATE AN ALIGNER
    aligner = Aligner(crs="EPSG:31370")
    # ADD A THEMATIC POLYGON ( a building-draft) TO THEMATIC DICTIONARY and LOAD into Aligner
    id = "my_building_id"
    thematic_dict = {
        id: geom_from_wkt(
            "POLYGON ((172355.21254837638116442 171905.62395133438985795, 172365.70286233740625903 171900.59734256137744524, 172359.91133483810699545 171887.48445011011790484, 172357.17948224407155067 171887.70299831763259135, 172355.86819299895432778 171884.42477520479587838, 172346.90771649059024639 171888.14009473266196437, 172355.21254837638116442 171905.62395133438985795))"
        )
    }

    loader = DictLoader(data_dict=thematic_dict)
    aligner.load_thematic_data(loader)

    building_tags = {"building": True}

    landuse_tags = {
        "landuse": [
            "residential",
            "commercial",
            "industrial",
            "retail",
            "farmland",
            "forest",
        ]
    }

    admin_tags = {"boundary": "administrative"}
    # Use a OSMLoader for the reference data
    loader = OSMLoader(osm_tags=landuse_tags, aligner=aligner)
    aligner.load_reference_data(loader)
    # EXECUTE THE ALIGNMENT
    relevant_distance = 5
    process_result = aligner.process(relevant_distance=relevant_distance)
    # PRINT RESULTS IN WKT
    print("result: " + process_result[id][relevant_distance]["result"].wkt)
    print(
        "added area: " + process_result[id][relevant_distance]["result_diff_plus"].wkt
    )
    print(
        "removed area: " + process_result[id][relevant_distance]["result_diff_min"].wkt
    )
