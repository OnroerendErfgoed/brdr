from datetime import date

import numpy as np

from brdr.aligner import Aligner
from brdr.enums import GRBType
from brdr.grb import GRBActualLoader
from brdr.grb import GRBFiscalParcelLoader
from brdr.grb import evaluate
from brdr.grb import get_geoms_affected_by_grb_change
from brdr.loader import DictLoader
from brdr.oe import OnroerendErfgoedLoader
from brdr.utils import get_series_geojson_dict

#from brdr.utils import merge_process_results

multi_as_single_modus = False

# Align the multipolygon to the fiscal parcels 2022

base_aligner = Aligner()
base_aligner.multi_as_single_modus = multi_as_single_modus
loader = OnroerendErfgoedLoader([9946])
base_aligner.load_thematic_data(loader)
base_year = "2022"
base_aligner.load_reference_data(
    GRBFiscalParcelLoader(year=base_year, aligner=base_aligner)
)
relevant_distance=2
base_process_result = base_aligner.process_dict_thematic(relevant_distance=relevant_distance)
#base_process_result = merge_process_results(base_process_result)
thematic_dict_formula = {}
thematic_dict_result = {}

# Create a dictionary with resulting geometries (aligned on Adpf2022) and a dictionary
# with the corresponding formula
for key in base_process_result:
    thematic_dict_result[key] = base_process_result[key][relevant_distance]["result"]
    thematic_dict_formula[key] = base_aligner.get_formula(thematic_dict_result[key])

# Determine all features that are possibly changed during timespan
base_aligner_result = Aligner()
base_aligner_result.multi_as_single_modus = multi_as_single_modus
base_aligner_result.load_thematic_data(DictLoader(thematic_dict_result))
dict_affected, dict_unchanged = get_geoms_affected_by_grb_change(
    base_aligner_result,
    grb_type=GRBType.ADP,
    date_start=date(2022, 1, 1),
    date_end=date.today(),
    one_by_one=False,
)
# Align the possibly affected geometry on the actual GRB parcels (evaluation)


actual_aligner = Aligner()
actual_aligner.multi_as_single_modus = multi_as_single_modus
loader = DictLoader(dict_affected)
actual_aligner.load_thematic_data(loader)
loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=actual_aligner)
actual_aligner.load_reference_data(loader)
series = np.arange(0, 200, 10, dtype=int) / 100
dict_series, dict_predicted, diffs_dict = actual_aligner.predictor(series)

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

print(fc["result"])
