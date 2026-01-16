from brdr.aligner import Aligner
from brdr.geometry_utils import geom_from_wkt
from brdr.loader import DictLoader

if __name__ == "__main__":
    """
    Example to load dat with a DictLoader, and also adding the properties to the result
    """
    # CREATE AN ALIGNER
    aligner = Aligner(crs="EPSG:31370")
    # ADD A THEMATIC POLYGON TO THEMATIC DICTIONARY and LOAD into Aligner
    id = "my_theme_id"
    thematic_dict = {id: geom_from_wkt("POLYGON ((0 0, 0 9, 5 10, 10 0, 0 0))")}
    # Add properties
    thematic_dict_properties = {
        id: {"propA": 1, "propB": 1.1, "propC": "dit is tekst", "propD": None}
    }
    loader = DictLoader(
        data_dict=thematic_dict, data_dict_properties=thematic_dict_properties
    )
    aligner.load_thematic_data(loader)
    # ADD A REFERENCE POLYGON TO REFERENCE DICTIONARY and LOAD into Aligner
    reference_dict = {"my_ref_id": geom_from_wkt("POLYGON ((0 1, 0 10,8 10,10 1,0 1))")}
    loader = DictLoader(reference_dict)
    aligner.load_reference_data(loader)
    # EXECUTE THE ALIGNMENT
    relevant_distance = 1
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
