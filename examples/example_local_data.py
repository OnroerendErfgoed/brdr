from brdr.aligner import Aligner
from brdr.enums import OpenbaarDomeinStrategy
from brdr.loader import GeoJsonFileLoader

if __name__ == "__main__":
    # Initiate brdr
    aligner = Aligner()
    # Load local thematic data and reference data
    loader = GeoJsonFileLoader(
        "../tests/testdata/themelayer_referenced.geojson", "id_theme"
    )
    aligner.load_thematic_data(loader)
    loader = GeoJsonFileLoader("../tests/testdata/reference_leuven.geojson", "capakey")
    aligner.load_reference_data(loader)
    # Example how to use the Aligner
    rel_dist = 1
    dict_results = aligner.process(
        relevant_distance=rel_dist,
        od_strategy=OpenbaarDomeinStrategy.SNAP_FULL_AREA_ALL_SIDE,
    )

    aligner.save_results("output/")
    # show_map(dict_results, aligner.dict_thematic, aligner.dict_reference)
