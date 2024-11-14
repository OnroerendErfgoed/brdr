from brdr.aligner import Aligner
from brdr.enums import AlignerInputType
from brdr.grb import update_to_actual_grb
from brdr.loader import GeoJsonFileLoader

aligner = Aligner()
# Load thematic data
loader = GeoJsonFileLoader("research.geojson", "OBJECTID")
aligner.load_thematic_data(loader)

fc = aligner.get_input_as_geojson(inputtype=AlignerInputType.THEMATIC)

fcs_actualisation = update_to_actual_grb(
    fc,
    id_theme_fieldname="OBJECTID",
    base_formula_field="brdr_formula",
    max_distance_for_actualisation=3,
    feedback=None,
)

print(fcs_actualisation)
