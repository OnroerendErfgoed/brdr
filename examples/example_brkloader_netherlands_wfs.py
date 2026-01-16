from brdr.aligner import Aligner
from brdr.geometry_utils import geom_from_wkt
from brdr.loader import DictLoader, WFSReferenceLoader
from brdr.nl.constants import BRK_CRS

if __name__ == "__main__":
    """
    Example to reference data with WFSReferenceLoader (Basisregistratie Kadaster - Netherlands): https://service.pdok.nl/kadaster/kadastralekaart/wfs/v5_0?request=GetCapabilities&service=WFS
    """
    # CREATE AN ALIGNER
    aligner = Aligner(crs=BRK_CRS)
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
    loader = WFSReferenceLoader(
        url="https://service.pdok.nl/kadaster/kadastralekaart/wfs/v5_0",
        id_property="identificatieLokaalID",
        typename="kadastralekaart:Perceel",
        aligner=aligner,
    )
    aligner.load_reference_data(loader)
    # EXECUTE THE ALIGNMENT
    relevant_distance = 5
    aligner_result = aligner.process(relevant_distances=[relevant_distance])
    # PRINT RESULTS IN WKT
    print("result: " + aligner_result.results[id][relevant_distance]["result"].wkt)
    print(
        "added area: "
        + aligner_result.results[id][relevant_distance]["result_diff_plus"].wkt
    )
    print(
        "removed area: "
        + aligner_result.results[id][relevant_distance]["result_diff_min"].wkt
    )
