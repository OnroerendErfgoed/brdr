from datetime import date

import numpy as np

from brdr.aligner import Aligner
from brdr.constants import EVALUATION_FIELD_NAME, FORMULA_FIELD_NAME
from brdr.enums import GRBType
from brdr.grb import (
    get_geoms_affected_by_grb_change,
    GRBFiscalParcelLoader,
    GRBActualLoader,
)
from brdr.loader import DictLoader
from brdr.oe import OnroerendErfgoedLoader
from brdr.utils import get_series_geojson_dict

base_aligner = Aligner()
loader = OnroerendErfgoedLoader([120288])
base_aligner.load_thematic_data(loader)
base_year = "2022"
base_aligner.load_reference_data(
    GRBFiscalParcelLoader(year=base_year, aligner=base_aligner)
)
relevant_distance = 3
base_process_result = base_aligner.process(relevant_distance=relevant_distance)
thematic_dict_formula = {}
thematic_dict_result = {}
for key in base_process_result:
    thematic_dict_result[key] = base_process_result[key][relevant_distance]["result"]
    thematic_dict_formula[key] = {
        FORMULA_FIELD_NAME: base_aligner.get_brdr_formula(thematic_dict_result[key])
    }
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

actual_aligner = Aligner()
actual_aligner.load_thematic_data(
    DictLoader(data_dict=dict_affected, data_dict_properties=thematic_dict_formula)
)
loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=actual_aligner)
actual_aligner.load_reference_data(loader)
series = np.arange(0, 300, 10, dtype=int) / 100

dict_evaluated, prop_dictionary = actual_aligner.compare(
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
for feature in fc["result"]["features"]:
    print(feature["properties"][EVALUATION_FIELD_NAME])
