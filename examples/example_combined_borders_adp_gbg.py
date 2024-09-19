from brdr.aligner import Aligner
from brdr.enums import GRBType
from brdr.grb import get_collection_grb_actual, GRBActualLoader
from brdr.loader import GeoJsonFileLoader, DictLoader
from brdr.utils import polygonize_reference_data, geojson_to_dicts
from examples import show_map, print_brdr_formula

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
        path_to_file="../tests/testdata/test_parcel_vs_building.geojson",
        id_property="theme_id",
    )
    aligner.load_thematic_data(loader)

    adploader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
    gbgloader = GRBActualLoader(grb_type=GRBType.GBG, partition=1000, aligner=aligner)
    collection_adp, name_reference_id_adp = get_collection_grb_actual(
        aligner.get_thematic_union(),
        grb_type=GRBType.ADP,
        partition=1000,
        date_start=None,
        date_end=None,
    )
    dict_adp, dict_adp_properties = geojson_to_dicts(
        collection_adp, name_reference_id_adp
    )

    collection_gbg, name_reference_id_gbg = get_collection_grb_actual(
        aligner.get_thematic_union(),
        grb_type=GRBType.GBG,
        partition=1000,
        date_start=None,
        date_end=None,
    )
    dict_gbg, dict_gbg_properties = geojson_to_dicts(
        collection_gbg, name_reference_id_gbg
    )

    dict_adp_gbg = dict_adp
    dict_adp_gbg.update(dict_gbg)  # combine 2 dictionaries
    # make a polygonized version of the reference data with non-overlapping polygons
    dict_ref = polygonize_reference_data(dict_adp_gbg)
    aligner.load_reference_data(DictLoader(dict_ref))

    rel_dist = 2
    dict_results = aligner.process(relevant_distances=[rel_dist], od_strategy=4)
    aligner.save_results("output/")
    show_map(dict_results, aligner.dict_thematic, aligner.dict_reference)
    print_brdr_formula(dict_results, aligner)
