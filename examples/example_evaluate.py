from datetime import date

import numpy as np
from shapely import from_wkt

from brdr.aligner import Aligner
from brdr.enums import GRBType
from brdr.grb import GRBActualLoader
from brdr.grb import GRBFiscalParcelLoader
from brdr.grb import evaluate
from brdr.grb import get_geoms_affected_by_grb_change
from brdr.loader import DictLoader, GeoJsonFileLoader
from brdr.utils import get_series_geojson_dict


def fid_to_geojson(geojson):
    fid = 1
    for f in geojson["features"]:
        f["properties"]["fid"] = str(fid)
        fid = fid + 1
        if f["geometry"]["type"] == "Polygon":
            f["geometry"] = {
                "type": "MultiPolygon",
                "coordinates": [f["geometry"]["coordinates"]],
            }

    return geojson


#
# thematic_dict = {
#     "theme_id_1": from_wkt(
#         "MultiPolygon (((174180.20077791667426936 171966.14649116666987538, "
# "174415.60530965600628406 171940.9636807945498731, "
# "174388.65236948925303295 171770.99678386366576888, "
# "174182.10876987033407204 171836.13745758961886168, "
# "174184.88916448061354458 171873.07698598300339654, "
# "174180.20077791667426936 171966.14649116666987538)))"
#     )
# }

# thematic_dict = {
#     "theme_id_1": from_wkt(
#         "MultiPolygon (((173463.11530961000244133 174423.83310307000647299, "
# "173460.22633100001257844 174422.02316300000529736, "
# "173455.24681099998997524 174429.98009100000490434, "
# "173454.4299790000077337 174429.34482699999352917, "
# "173452.06690700000035577 174432.43058700000983663, "
# "173451.25743500000680797 174431.8672589999914635, "
# "173448.74844299998949282 174434.96249100001296028, "
# "173448.5809550000121817 174435.80485899999621324, "
# "173454.82841871041455306 174442.46780387416947633, "
# "173461.44169100001454353 174446.50898700000834651, "
# "173472.15932299999985844 174429.49919500001124106, "
# "173466.18524341000011191 174425.75641125999391079, "
# "173466.9701960513193626 174424.8217541387421079, "
# "173462.59915620859828778 174424.8217541387421079, "
# "173463.11530961000244133 174423.83310307000647299)))"
#     )
# }
# Polygon ((173455.24681099998997524 174429.9801549999974668, "
# "173454.4299790000077337 174429.34482699999352917, "
# "173452.06690700000035577 174432.43058700000983663, "
# "173451.25743500000680797 174431.8672589999914635, "
# "173448.74844299998949282 174434.96249100001296028, "
# "173448.5809550000121817 174435.80485899999621324, "
# "173455.39772300000186078 174441.47852299999794923, "
# "173461.44169100001454353 174446.50898700000834651, "
# "173472.15932299999985844 174429.49919500001124106, "
# "173466.18524300001445226 174425.75641100000939332, "
# "173460.22633100001257844 174422.02316300000529736, "
# "173455.24681099998997524 174429.9801549999974668))"

# MultiPolygon (((173463.11530961000244133 174423.83310307000647299, "
# "173460.22633100001257844 174422.02316300000529736, "
# "173455.24681099998997524 174429.98009100000490434, "
# "173454.4299790000077337 174429.34482699999352917, "
# "173452.06690700000035577 174432.43058700000983663, "
# "173451.25743500000680797 174431.8672589999914635, "
# "173448.74844299998949282 174434.96249100001296028, "
# "173448.5809550000121817 174435.80485899999621324, "
# "173454.82841871041455306 174442.46780387416947633, "
# "173461.44169100001454353 174446.50898700000834651, "
# "173472.15932299999985844 174429.49919500001124106, "
# "173466.18524341000011191 174425.75641125999391079, "
# "173466.9701960513193626 174424.8217541387421079, "
# "173462.59915620859828778 174424.8217541387421079, "
# "173463.11530961000244133 174423.83310307000647299)))"

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
base_process_result = base_aligner.process_dict_thematic(
    relevant_distance=relevant_distance
)
thematic_dict_formula = {}
thematic_dict_result = {}
for key in base_process_result:
    thematic_dict_result[key] = base_process_result[key][relevant_distance]["result"]
    thematic_dict_formula[key] = base_aligner.get_formula(thematic_dict_result[key])
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
actual_aligner.load_thematic_data(DictLoader(dict_affected))
actual_aligner.load_reference_data(
    GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=actual_aligner)
)
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
print(fc["result"])
fcs = actual_aligner.get_series_as_geojson(formula=True)
print(fcs["result"])

for feature in fc["result"]["features"]:
    print(
        feature["properties"][actual_aligner.name_thematic_id]
        + ": "
        + feature["properties"]["evaluation"]
    )

geojson = fid_to_geojson(fc["result"])


print(geojson)
