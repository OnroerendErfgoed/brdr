from brdr.aligner import Aligner
from examples import plot_series
from examples import show_results

if __name__ == "__main__":
    # Initiate brdr
    aligner = Aligner()
    # Load thematic data
    # x.load_thematic_data_file("../tests/testdata/theme_leuven.geojson", 'aanduid_id')
    # x.load_thematic_data_file("../tests/testdata/theme.geojson", 'theme_identifier')
    aligner.load_thematic_data_file(
        "../tests/testdata/themelayer_referenced.geojson", "id_theme"
    )
    # x.load_thematic_data_file("../tests/testdata/zwin.geojson", 'id')
    # x.load_thematic_data_file("theme_vmm.geojson", 'unique_id')

    # gebruik de actuele adp-percelen adp= administratieve percelen
    aligner.load_reference_data_grb_actual(grb_type="adp", partition=1000)

    # # gebruik de GRB-gebouwen, gbg= gebouw aan de grond
    # x.load_reference_data_grb_actual('gbg')

    # x.load_reference_data_file("../tests/testdata/reference_leuven.geojson", 'capakey')
    # x.load_reference_data_file("reference_vmm.geojson", 'ref_identifier')

    # Example how to use the Aligner
    r, rd, rd_plus, rd_min, sd, si = aligner.process_dict_thematic(10, 3)
    aligner.export_results("output/")
    show_results(r, rd_plus,rd_min,aligner.dict_thematic, aligner.dict_reference)

    r, rd, rd_plus, rd_min, sd, si = aligner.process_dict_thematic(6, 2)
    aligner.export_results("output/")
    show_results(r, rd_plus,rd_min,aligner.dict_thematic, aligner.dict_reference)
    # for key in r:
    #     x.get_formula(r[key])

    # Example how to use a series (for histogram)
    # series = [0.2, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 6, 8, 10]
    series = [0.1, 0.2, 0.3, 0.4, 0.5, 1, 2]
    resulting_areas = aligner.process_series(series, 2, 50)
    plot_series(series, resulting_areas)

    # Example how to use the Aligner with full_overlap_percentage=-1 (original
    # border will be used)
    r, rd, rd_plus, rd_min, sd, si = aligner.process_dict_thematic(30, 1, -1)
    aligner.export_results("output/")
    show_results(r, rd_plus,rd_min,aligner.dict_thematic, aligner.dict_reference)
    # for key in r:
    #     x.get_formula(r[key])