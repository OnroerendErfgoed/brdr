from brdr.aligner import Aligner
from brdr.enums import OpenbaarDomeinStrategy
from examples import show_map

if __name__ == "__main__":
    # Initiate brdr
    aligner = Aligner()
    #Load local thematic data and reference data
    aligner.load_thematic_data_file(
        "../tests/testdata/themelayer_referenced.geojson", "id_theme"
    )
    aligner.load_reference_data_file('../tests/testdata/reference_leuven.geojson', "capakey")
    # Example how to use the Aligner
    rel_dist = 1
    dict_results_by_distance = {}
    dict_results_by_distance[aligner.relevant_distance] = aligner.process_dict_thematic(relevant_distance=rel_dist, od_strategy=OpenbaarDomeinStrategy.SNAP_FULL_AREA_ALL_SIDE)
    aligner.export_results("output/")
    #show_map(dict_results_by_distance, aligner.dict_thematic, aligner.dict_reference)
