from brdr.enums import GRBType
from brdr.grb import get_collection_grb_actual
from brdr.loader import GRBActualLoader, GeoJsonFileLoader
from brdr.aligner import Aligner
from brdr.utils import polygonize_reference_data, collection_to_dict
from examples import show_map, print_formula

# example to test what happens if we combine borders
# (so thematic data can use both polygons)

# If we just combine sets of polygons (fe parcels and buildings) these polygons will
# overlap, and gives some unwanted/unexpected results. The reference data can be
# preprocessed with 'shapely-polygonize', to create a new set of non-overlapping
# polygons. These non-overlapping polygons can be used as reference-data to align the
# theme-data.


if __name__ == "__main__":
    # Initiate brdr
    aligner = Aligner()

    # Load thematic data & reference data
    loader = GeoJsonFileLoader(
        path_to_file="../tests/testdata/test_parcel_vs_building.geojson", id_property="theme_id"
    )
    aligner.load_thematic_data(loader)

    adploader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
    gbgloader = GRBActualLoader(grb_type=GRBType.GBG, partition=1000, aligner=aligner)
    collection_adp, name_reference_id_adp = get_collection_grb_actual(aligner._get_thematic_union(), grb_type=GRBType.ADP, partition=1000,
                              date_start=None, date_end=None)
    dict_adp = collection_to_dict(collection_adp,name_reference_id_adp)

    collection_gbg, name_reference_id_gbg = get_collection_grb_actual(aligner._get_thematic_union(), grb_type=GRBType.GBG, partition=1000,
                              date_start=None, date_end=None)
    dict_gbg = collection_to_dict(collection_gbg,name_reference_id_gbg)

    dict_adp_gbg = dict_adp
    dict_adp_gbg.update(dict_gbg)  # combine 2 dictionaries
    # make a polygonized version of the reference data with non-overlapping polygons
    dict_ref = polygonize_reference_data(dict_adp_gbg)
    aligner.load_reference_data_dict(dict_ref)

    rel_dist = 2
    dict_results_by_distance = {}
    dict_results_by_distance[rel_dist] = aligner.process_dict_thematic(rel_dist, 4)
    aligner.export_results("output/")
    show_map(dict_results_by_distance, aligner.dict_thematic, aligner.dict_reference)
    print_formula(dict_results_by_distance, aligner)
