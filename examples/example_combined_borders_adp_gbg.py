from brdr.aligner import Aligner
from brdr.be.grb.enums import GRBType
from brdr.be.grb.loader import GRBActualLoader
from brdr.enums import AlignerResultType
from brdr.loader import GeoJsonFileLoader, DictLoader
from brdr.utils import polygonize_reference_data
from brdr.viz import show_map, print_observation_of_aligner_results

if __name__ == "__main__":
    """
    # example to test what happens if we combine borders
    # (so thematic data can use both polygons)

    # If we just combine sets of polygons (fe parcels and buildings) these polygons will
    # overlap, and gives some unwanted/unexpected results. The reference data can be
    # preprocessed with 'shapely-polygonize', to create a new set of non-overlapping
    # polygons. These non-overlapping polygons can be used as reference-data to align the
    # theme-data.
    """
    # Initiate brdr
    aligner = Aligner()

    # Load thematic data & reference data
    loader = GeoJsonFileLoader(
        path_to_file="../tests/testdata/test_parcel_vs_building.geojson",
        id_property="theme_id",
    )
    aligner.load_thematic_data(loader)

    adp_loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
    gbg_loader = GRBActualLoader(grb_type=GRBType.GBG, partition=1000, aligner=aligner)
    reference_data_adp = adp_loader.load_data()
    reference_data_gbg = gbg_loader.load_data()

    adp_features = reference_data_adp.features
    gbg_features = reference_data_adp.features
    # Create new ref dictionary based on 2 reference sources
    dict_ref = {}
    for ref_id, feature in adp_features.items():
        dict_ref[ref_id] = feature.geometry
    for ref_id, feature in gbg_features.items():
        dict_ref[ref_id] = feature.geometry
    # make a polygonized version of the reference data with non-overlapping polygons
    dict_ref = polygonize_reference_data(dict_ref)
    aligner.load_reference_data(DictLoader(dict_ref))

    rel_dist = 2
    aligner_result = aligner.process(relevant_distances=[rel_dist])
    aligner_result.save_results(path="output/", aligner=aligner)
    thematic_geometries = {
        key: feat.geometry for key, feat in aligner.thematic_data.features.items()
    }
    reference_geometries = {
        key: feat.geometry for key, feat in aligner.reference_data.features.items()
    }
    show_map(aligner_result.results, thematic_geometries, reference_geometries)
    print_observation_of_aligner_results(aligner_result.results, aligner)

    fcs = aligner_result.get_results_as_geojson(
        result_type=AlignerResultType.PROCESSRESULTS, aligner=aligner
    )
    print(fcs["result"])

    print(aligner.thematic_data.to_geojson())
