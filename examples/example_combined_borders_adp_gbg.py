from brdr.aligner import Aligner
from examples import show_results

# example to test what happens if we combine borders
# (so thematic data can use both polygons)

# If we just combine sets of polygons (fe parcels and buildings) these polygons will
# overlap, and gives some unwanted/unexpected results. The reference data can be
# preprocessed with 'shapely-polygonize', to create a new set of non-overlapping
# polygons. These non-overlapping polygons can be used as reference-data to align the
# theme-data.


if __name__ == "__main__":
    # Initiate brdr
    x = Aligner()

    # Load thematic data & reference data
    x.load_thematic_data_file(
        "../tests/testdata/test_parcel_vs_building.geojson", "theme_id"
    )
    dict_adp, name_reference_id_adp = x.get_reference_data_dict_grb_actual(
        grb_type="adp", partition=1000
    )
    dict_gbg, name_reference_id_gbg = x.get_reference_data_dict_grb_actual(
        grb_type="gbg", partition=1000
    )
    dict_adp_gbg = dict_adp
    dict_adp_gbg.update(dict_gbg)  # combine 2 dictionaries
    # make a polygonized version of the reference data with non-overlapping polygons
    dict_ref = dict_adp_gbg
    x.load_reference_data_dict(dict_ref)

    r, rd, rd_plus, rd_min, sd, si = x.process_dict_thematic(2, 4)
    x.export_results("output/")
    show_results(r, rd)

    for key in r:
        x.get_formula(r[key])
