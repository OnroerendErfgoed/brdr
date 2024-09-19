from datetime import date

import numpy as np
from shapely import from_wkt

from brdr.aligner import Aligner
from brdr.constants import FORMULA_FIELD_NAME, EVALUATION_FIELD_NAME
from brdr.enums import GRBType
from brdr.grb import GRBActualLoader
from brdr.grb import GRBFiscalParcelLoader
from brdr.grb import get_geoms_affected_by_grb_change
from brdr.loader import DictLoader, GeoJsonFileLoader
from brdr.utils import get_series_geojson_dict

thematic_dict = {
    "theme_id_1": from_wkt(
        "Polygon ((174072.91453437806922011 179188.47430499014444649, 174121.17416846146807075 179179.98909460185677744, 174116.93156326730968431 179156.47799081765697338, 174110.56765547610120848 179152.58893605635967106, 174069.37903004963300191 179159.30639428040012717, 174069.37903004963300191 179159.30639428040012717, 174070.97000699743512087 179169.7361320493509993, 174072.91453437806922011 179188.47430499014444649))"
    )
}
base_aligner = Aligner()
loader = GeoJsonFileLoader("themelayer.geojson", "theme_identifier")
base_aligner.load_thematic_data(loader)
base_year = "2022"
base_aligner.load_reference_data(
    GRBFiscalParcelLoader(year=base_year, aligner=base_aligner)
)
relevant_distance = 2
base_process_result = base_aligner.process(relevant_distance=relevant_distance)
thematic_dict_formula = {}
thematic_dict_result = {}
for key in base_process_result:
    thematic_dict_result[key] = base_process_result[key][relevant_distance]["result"]
    thematic_dict_formula[key] = {
        FORMULA_FIELD_NAME: base_aligner.get_brdr_formula(thematic_dict_result[key])
    }
    print(key + ": " + thematic_dict_result[key].wkt)
    print(key + ": " + str(thematic_dict_formula[key]))
base_aligner_result = Aligner()
base_aligner_result.load_thematic_data(DictLoader(thematic_dict_result))
dict_affected, dict_unchanged = get_geoms_affected_by_grb_change(
    base_aligner_result,
    grb_type=GRBType.ADP,
    date_start=date(2022, 1, 1),
    date_end=date.today(),
    one_by_one=False,
)
if dict_affected == {}:
    print("No affected dicts")
    exit()
for key, value in dict_affected.items():
    print(key + ": " + value.wkt)
actual_aligner = Aligner()
loader = DictLoader(dict_affected)
actual_aligner.load_thematic_data(
    DictLoader(data_dict=dict_affected, data_dict_properties=thematic_dict_formula)
)
actual_aligner.load_reference_data(
    GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=actual_aligner)
)
actual_aligner.relevant_distances = np.arange(0, 200, 10, dtype=int) / 100
dict_evaluated, prop_dictionary = actual_aligner.compare(
    # thematic_dict_formula=thematic_dict_formula,
    threshold_area=5,
    threshold_percentage=1,
    dict_unchanged=dict_unchanged,
)

fc = get_series_geojson_dict(
    dict_evaluated,
    crs=actual_aligner.CRS,
    id_field=actual_aligner.name_thematic_id,
    series_prop_dict=prop_dictionary,
)
print(fc["result"])
fcs = actual_aligner.get_results_as_geojson(formula=True)
print(fcs["result"])

for feature in fc["result"]["features"]:
    print(
        feature["properties"][actual_aligner.name_thematic_id]
        + ": "
        + feature["properties"][EVALUATION_FIELD_NAME]
    )
