from brdr.aligner import Aligner
from brdr.geometry_utils import geom_from_wkt
from brdr.loader import DictLoader

THEME_ID = "theme_id_1"
REFERENCE_ID = "ref_id_1"
RELEVANT_DISTANCE = 1

aligner = Aligner(crs="EPSG:31370")

# Load thematic geometry.
thematic_geometries = {
    THEME_ID: geom_from_wkt("POLYGON ((0 0, 0 9, 5 10, 10 0, 0 0))")
}
aligner.load_thematic_data(DictLoader(thematic_geometries))

# Load reference geometry.
reference_geometries = {
    REFERENCE_ID: geom_from_wkt("POLYGON ((0 1, 0 10,8 10,10 1,0 1))")
}
aligner.load_reference_data(DictLoader(reference_geometries))

# Run alignment for a single relevant distance.
aligner_result = aligner.process(relevant_distances=[RELEVANT_DISTANCE])
process_results = aligner_result.get_results(aligner=aligner)

# Print result geometries in WKT.
print("result: " + process_results[THEME_ID][RELEVANT_DISTANCE]["result"].wkt)
print(
    "added area: "
    + process_results[THEME_ID][RELEVANT_DISTANCE]["result_diff_plus"].wkt
)
print(
    "removed area: "
    + process_results[THEME_ID][RELEVANT_DISTANCE]["result_diff_min"].wkt
)
