from datetime import date

import numpy as np

from brdr.aligner import Aligner
from brdr.constants import EVALUATION_FIELD_NAME, FORMULA_FIELD_NAME
from brdr.enums import GRBType, AlignerResultType
from brdr.grb import (
    get_affected_by_grb_change,
    GRBFiscalParcelLoader,
    GRBActualLoader,
)
from brdr.loader import DictLoader
from brdr.oe import OnroerendErfgoedLoader

base_aligner = Aligner()
loader = OnroerendErfgoedLoader([120288,120108])
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
    print(key + ": " + thematic_dict_result[key].wkt)
    print(key + ": " + str(thematic_dict_formula[key]))
base_aligner_result = Aligner()
base_aligner_result.load_thematic_data(DictLoader(thematic_dict_result))
affected, unaffected = get_affected_by_grb_change(
    dict_thematic = thematic_dict_result,
    grb_type=GRBType.ADP,
    date_start=date(2022, 1, 1),
    date_end=date.today(),
    one_by_one=False,
)
if len(affected)==0:
    print("No affected dicts")
    exit()
print("Affected_IDs: " + str(affected))
actual_aligner = Aligner()
actual_aligner.load_thematic_data(
    DictLoader(data_dict=thematic_dict_result, data_dict_properties=thematic_dict_formula)
)
actual_aligner.load_reference_data(
    GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=actual_aligner)
)
actual_aligner.relevant_distances = np.arange(0, 200, 10, dtype=int) / 100
dict_evaluated, prop_dictionary = actual_aligner.compare(ids_to_compare=affected)

fc = actual_aligner.get_results_as_geojson(resulttype=AlignerResultType.EVALUATED_PREDICTIONS)
print(fc["result"])

for feature in fc["result"]["features"]:
    id = feature["properties"][actual_aligner.name_thematic_id]
    evaluation = feature["properties"][EVALUATION_FIELD_NAME]
    print(id + ": " + evaluation)