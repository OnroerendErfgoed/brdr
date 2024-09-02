from datetime import date

import numpy as np

from brdr.aligner import Aligner
from brdr.enums import GRBType
from brdr.grb import (
    get_geoms_affected_by_grb_change,
    evaluate,
    GRBFiscalParcelLoader,
    GRBActualLoader,
)
from brdr.loader import DictLoader
from brdr.utils import get_series_geojson_dict, get_oe_dict_by_ids

# dict_theme = get_oe_dict_by_ids([125610,148305,127615,122316,120153,124699,115489,
# 120288,120387,124762,148143,116141])
dict_theme = get_oe_dict_by_ids([10047, 10048, 10049, 10050, 10051, 10056])
print(dict_theme)

base_aligner = Aligner()
base_aligner.load_thematic_data(DictLoader(dict_theme))
base_year = "2022"
base_aligner.load_reference_data(
    GRBFiscalParcelLoader(year=base_year, aligner=base_aligner)
)
base_process_result = base_aligner.process_dict_thematic(relevant_distance=2)
thematic_dict_formula = {}
thematic_dict_result = {}
for key in base_process_result:
    thematic_dict_result[key] = base_process_result[key]["result"]
    thematic_dict_formula[key] = base_aligner.get_formula(thematic_dict_result[key])
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
loader = DictLoader(dict_affected)
actual_aligner.load_thematic_data(loader)
loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=actual_aligner)
actual_aligner.load_reference_data(loader)
series = np.arange(0, 200, 10, dtype=int) / 100
dict_series, dict_predicted, diffs_dict = actual_aligner.predictor(series)

# diffs_dict=merge_diffs_dict(diffs_dict)

dict_evaluated, prop_dictionary = evaluate(
    actual_aligner,
    dict_series,
    dict_predicted,
    thematic_dict_formula,
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
    print(feature["properties"]["evaluation"])
