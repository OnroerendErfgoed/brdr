from brdr.aligner import Aligner
from brdr.be.be import BeCadastralParcelLoader
from brdr.geometry_utils import geom_from_wkt
from brdr.loader import DictLoader

if __name__ == "__main__":
    """
    Example to reference data with CadastralLoader (FODFIN Cadastral Parcels - Belgium)
    """
    aligner = Aligner(crs="EPSG:3812")
    # ADD A THEMATIC POLYGON ( a building-draft) TO THEMATIC DICTIONARY and LOAD into Aligner
    id = "my_countour_id"
    thematic_dict = {
        id: geom_from_wkt(
            "POLYGON ((672875.34840881277341396 671536.28253797488287091, 672920.61456579703371972 671524.15767449699342251, 672897.17316307302098721 671459.49173594801686704, 672859.7208069966873154 671472.15548224712256342, 672875.34840881277341396 671536.28253797488287091))"
        )
    }

    loader = DictLoader(data_dict=thematic_dict)
    aligner.load_thematic_data(loader)

    # Use BeCadastralParcel for the reference data
    loader = BeCadastralParcelLoader(aligner=aligner)
    aligner.load_reference_data(loader)
    # EXECUTE THE ALIGNMENT
    relevant_distance = 5
    aligner_result = aligner.process(relevant_distances=[relevant_distance])
    process_results = aligner_result.get_results(aligner=aligner)
    # PRINT RESULTS IN WKT
    print("result: " + process_results[id][relevant_distance]["result"].wkt)
    print(
        "added area: " + process_results[id][relevant_distance]["result_diff_plus"].wkt
    )
    print(
        "removed area: " + process_results[id][relevant_distance]["result_diff_min"].wkt
    )
