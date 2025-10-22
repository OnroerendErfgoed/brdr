from brdr.aligner import Aligner
from brdr.geometry_utils import geom_from_wkt
from brdr.loader import DictLoader
from brdr.loaders.nl import BRKLoader, BRKType
from brdr.osm import OSMLoader

if __name__ == "__main__":
    """
    Example to reference data with BRKLoader (Basisregistratie Kadaster - Netherlands)
    """
    # CREATE AN ALIGNER
    aligner = Aligner(crs="EPSG:28992")
    # ADD A THEMATIC POLYGON ( a building-draft) TO THEMATIC DICTIONARY and LOAD into Aligner
    id = "my_countour_id"
    thematic_dict = {
        id: geom_from_wkt(
            "POLYGON ((113607.2723651389387669 404726.03122497239382938, 113733.67886830697534606 404725.35525436722673476, 113723.53930922932340764 404532.70363189186900854, 113636.33910116153128911 404540.81527915399055928, 113607.2723651389387669 404726.03122497239382938))"
                )
    }

    loader = DictLoader(data_dict=thematic_dict)
    aligner.load_thematic_data(loader)

    # Use BRKLoader for the reference data
    loader = BRKLoader(brk_type=BRKType.perceel, partition=1000, aligner=aligner)
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
