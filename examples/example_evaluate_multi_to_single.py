from datetime import date

import numpy as np

from brdr.aligner import Aligner
from brdr.enums import GRBType
from brdr.grb import GRBActualLoader
from brdr.grb import GRBFiscalParcelLoader
from brdr.grb import evaluate
from brdr.grb import get_geoms_affected_by_grb_change
from brdr.loader import DictLoader
from brdr.utils import get_oe_dict_by_ids
from brdr.utils import get_series_geojson_dict
from brdr.utils import merge_process_results
from brdr.utils import multipolygons_to_singles

# thematic_dict = {
#     1: from_wkt(
#         "MultiPolygon (((173503.67240000001038425 179630.29370000000926666, "
#         "173517.47036350998678245 179636.5546227699960582, "
#         "173526.74985151999862865 179640.51110276998952031, "
#         "173534.37954751998768188 179620.73651076000533067, "
#         "173538.27292352000949904 179622.00793476001126692, "
#         "173540.02590000000782311 179616.82500000001164153, "
#         "173540.5209000000031665 179615.20139999999082647, "
#         "173541.88829999999143183 179610.71599999998579733, "
#         "173545.56226753001101315 179598.47103874001186341, "
#         "173535.27874752000207081 179595.28806274000089616, "
#         "173524.47042751000844873 179591.97017474001040682, "
#         "173507.03150000001187436 179623.90200000000186265, "
#         "173504.75964350000140257 179629.70476677000988275, "
#         "173503.67240000001038425 179630.29370000000926666)),"
#         "((173606.28779999999096617 179670.28599999999278225, "
#         "173610.20720000000437722 179671.40659999998752028, "
#         "173614.6045999999914784 179672.66380000000935979, "
#         "173631.59880000000703149 179677.52259999999660067, "
#         "173631.89120000001275912 179677.60620000000926666, "
#         "173633.32480000000214204 179678.01610000000800937, "
#         "173633.32519999999203719 179678.01610000000800937, "
#         "173636.02556358999572694 179678.24044679998769425, "
#         "173637.36501959001179785 179678.00627079998957925, "
#         "173638.8345875900122337 179677.26476680001360364, "
#         "173639.96585159000824206 179676.19820680000702851, "
#         "173640.61596359001123346 179675.11871879998943768, "
#         "173641.02953159000026062 179673.73823879999690689, "
#         "173641.81839559000218287 179668.88179079000838101, "
#         "173643.26172358999610879 179659.99692678998690099, "
#         "173644.17860000001383014 179654.35260000001289882, "
#         "173644.70509999999194406 179651.11160000000381842, "
#         "173646.01050000000395812 179643.0746999999973923, "
#         "173646.98720000000321306 179637.06119999999646097, "
#         "173647.33480000001145527 179634.92029999999795109, "
#         "173648.38750000001164153 179628.81059999999706633, "
#         "173648.38829999999143183 179628.80609999998705462, "
#         "173648.63229999999748543 179627.39060000001336448, "
#         "173648.78820000000996515 179626.48589999999967404, "
#         "173637.05239999998593703 179624.23920000001089647, "
#         "173628.10649999999441206 179622.52669999998761341, "
#         "173626.2447999999858439 179622.17029999999795109, "
#         "173623.44330000001355074 179621.6339999999909196, "
#         "173612.57180000000516884 179619.55280000000493601, "
#         "173609.84570000000530854 179627.77809999999590218, "
#         "173606.07050000000162981 179639.75099999998928979, "
#         "173602.54409999999916181 179650.93530000001192093, "
#         "173600.03229999999166466 179658.901400000002468, "
#         "173599.82219999999506399 179659.60740000000805594, "
#         "173598.40669999999227002 179664.36720000000786968, "
#         "173597.43900000001303852 179667.62049999998998828, "
#         "173597.43859999999403954 179667.62040000001434237, "
#         "173597.401400000002468 179667.74530000000959262, "
#         "173606.28779999999096617 179670.28599999999278225)))"
#     )
# }
thematic_dict = get_oe_dict_by_ids([9946])
# Align the multipolygon to the fiscal parcels 2022
# thematic_dict = multipolygons_to_singles(thematic_dict)
base_aligner = Aligner()
base_aligner.multi_as_single_modus = False
base_aligner.load_thematic_data(DictLoader(thematic_dict))
base_year = "2022"
base_aligner.load_reference_data(
    GRBFiscalParcelLoader(year=base_year, aligner=base_aligner)
)
base_process_result = base_aligner.process_dict_thematic(relevant_distance=2)
base_process_result = merge_process_results(base_process_result)
thematic_dict_formula = {}
thematic_dict_result = {}

# Create a dictionary with resulting geometries (aligned on Adpf2022) and a dictionary
# with the corresponding formula
for key in base_process_result:
    thematic_dict_result[key] = base_process_result[key]["result"]
    thematic_dict_formula[key] = base_aligner.get_formula(thematic_dict_result[key])

thematic_dict_result = multipolygons_to_singles(thematic_dict_result)
# Determine all features that are possibly changed during timespan
base_aligner_result = Aligner()
# base_aligner.multi_as_single_modus=False
base_aligner_result.load_thematic_data(DictLoader(thematic_dict_result))
dict_affected, dict_unchanged = get_geoms_affected_by_grb_change(
    base_aligner_result,
    grb_type=GRBType.ADP,
    date_start=date(2022, 1, 1),
    date_end=date.today(),
    one_by_one=False,
)
# Align the possibly affected geometry on the actual GRB parcels (evaluation)

# dict_affected = multipolygons_to_singles(dict_affected)
actual_aligner = Aligner()
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
