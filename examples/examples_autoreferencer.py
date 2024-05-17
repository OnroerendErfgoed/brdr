from brdr.auto_referencer import AutoReferencer
from examples import plot_diffs
from examples import show_results

if __name__ == "__main__":
    # Initiate brdr
    auto_referencer = AutoReferencer()
    # Load thematic data
    # x.load_thematic_data_file("../tests/testdata/theme_leuven.geojson", 'aanduid_id')
    # x.load_thematic_data_file("../tests/testdata/theme.geojson", 'theme_identifier')
    auto_referencer.load_thematic_data_file(
        "../tests/testdata/themelayer_referenced.geojson", "id_theme"
    )
    # x.load_thematic_data_file("../tests/testdata/zwin.geojson", 'id')
    # x.load_thematic_data_file("theme_vmm.geojson", 'unique_id')

    # gebruik de actuele adp-percelen adp= administratieve percelen
    auto_referencer.load_reference_data_grb_actual(grb_type="adp", partition=1000)

    # # gebruik de GRB-gebouwen, gbg= gebouw aan de grond
    # x.load_reference_data_grb_actual('gbg')

    # x.load_reference_data_file("../tests/testdata/reference_leuven.geojson", 'capakey')
    # x.load_reference_data_file("reference_vmm.geojson", 'ref_identifier')

    # Example how to use the AutoReferencer
    r, rd, *_ = auto_referencer.process_dict_thematic(10, 3)
    auto_referencer.export_results("output/")
    show_results(r, rd)

    r, rd, *_ = auto_referencer.process_dict_thematic(6, 2)
    auto_referencer.export_results("output/")
    show_results(r, rd)
    # for key in r:
    #     x.get_formula(r[key])

    # Example how to use a series (for histogram)
    # series = [0.2, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 6, 8, 10]
    series = [0.1, 0.2, 0.3, 0.4, 0.5, 1, 2]
    resulting_areas = auto_referencer.process_series(series, 2, 50)
    plot_diffs(series, resulting_areas)

    # Example how to use the AutoReferencer with full_overlap_percentage=-1 (original
    # border will be used)
    r, rd, rd_plus, rd_min, sd, si = auto_referencer.process_dict_thematic(30, 1, -1)
    auto_referencer.export_results("output/")
    show_results(r, rd)
    # for key in r:
    #     x.get_formula(r[key])

    # #Example how to use the ThemeFeature
    # p = shapely.wkt.loads('MultiPolygon(
    # ((173958.08938021279755048 174805.31052483833627775,
    # 173965.8578134463459719 174791.03602877166122198,
    # 173977.60756871209014207 174798.22182951270951889,
    # 173982.46283948307973333 174788.89970963244559243,
    # 173990.32837813204969279 174792.97813708006287925,
    # 173977.41335788124706596 174815.99212053447263315,
    # 173958.08938021279755048 174805.31052483833627775)))'
    # )
    # p = shapely.wkt.loads(
    #     "MultiPolygon("
    #     "(("
    #     "174523.92973850725684315 170552.32748373309732415, "
    #     "174665.33266440019360743 170520.55154533020686358, "
    #     "174666.92146132033667527 170460.17726236468297429, "
    #     "174525.51853542739991099 170490.36440384743036702, "
    #     "174523.92973850725684315 170552.32748373309732415)))"
    # )
    # tf = ThemeFeature(p, "id1")
    # tf.set_auto_referencer(x)
    # tf.align(2, -1, 50)
    # print(tf.geometry_result)
