
from datetime import date

import numpy as np
from shapely import from_wkt
from brdr.aligner import Aligner
from brdr.enums import GRBType
from brdr.geometry_utils import get_bbox
from brdr.grb import (
    get_geoms_affected_by_grb_change, evaluate, get_collection_grb_fiscal_parcels,
)
from brdr.loader import DictLoader, GeoJsonLoader, GRBFiscalParcelLoader
from brdr.loader import GRBActualLoader
from brdr.utils import get_series_geojson_dict

thematic_dict = {
    "theme_id_1": from_wkt(
        "MultiPolygon (((174180.20077791667426936 171966.14649116666987538, 174415.60530965600628406 171940.9636807945498731, 174388.65236948925303295 171770.99678386366576888, 174182.10876987033407204 171836.13745758961886168, 174184.88916448061354458 171873.07698598300339654, 174180.20077791667426936 171966.14649116666987538)))"
    )
}
bbox = get_bbox(thematic_dict["theme_id_1"])
base_aligner = Aligner()
base_aligner.load_thematic_data(DictLoader(thematic_dict))
base_year ="2022"
#collection_fiscal_parcels = get_collection_grb_fiscal_parcels(base_year, bbox=bbox)
#base_aligner.load_reference_data(GeoJsonLoader(collection_fiscal_parcels, "CAPAKEY"))
base_aligner.load_reference_data(GRBFiscalParcelLoader(year=base_year,aligner =base_aligner))
base_process_result = base_aligner.process_dict_thematic(relevant_distance=1)
thematic_dict_formula = {}
for key in base_process_result:
    thematic_dict[key] = base_process_result[key]["result"]
    thematic_dict_formula[key] = base_aligner.get_formula(thematic_dict[key])
dict_affected = get_geoms_affected_by_grb_change(
    thematic_dict,
    grb_type=GRBType.ADP,
    date_start=date(2022, 1, 1),
    date_end=date.today(),
    one_by_one=False,
)

series = np.arange(0, 200, 10, dtype=int) / 100

actual_aligner = Aligner()
loader = DictLoader(dict_affected)
actual_aligner.load_thematic_data(loader)
loader = GRBActualLoader(grb_type=GRBType.ADP, partition=0, aligner=actual_aligner)
actual_aligner.load_reference_data(loader)



dict_evaluated,prop_dictionary = evaluate(actual_aligner,thematic_dict_formula, series, )
fc = get_series_geojson_dict(
            dict_evaluated,
            crs=actual_aligner.CRS,
            id_field=actual_aligner.name_thematic_id,
            series_prop_dict=prop_dictionary,
        )

print (fc)