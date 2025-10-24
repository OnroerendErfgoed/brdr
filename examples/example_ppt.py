from brdr.aligner import Aligner
from brdr.be.grb.enums import GRBType
from brdr.be.grb.loader import GRBActualLoader
from brdr.enums import AlignerResultType
from brdr.geometry_utils import geom_from_wkt
from brdr.loader import DictLoader

# CREATE AN ALIGNER
aligner = Aligner(crs="EPSG:31370")
# Load thematic WKT into Aligner
thematic_data = {
    "theme_id_1": geom_from_wkt("POLYGON ((170996.19647056088433601 170611.47279253022861667, 171036.0500100600766018 170583.59100930689601228, 171025.62942813205881976 170566.12672797418781556, 170984.43460815289290622 170595.98966291974647902, 170996.19647056088433601 170611.47279253022861667))")
}
thematic_loader = DictLoader(thematic_data)
aligner.load_thematic_data(thematic_loader)

# Load reference data (GRB parcels)
reference_loader = GRBActualLoader(grb_type=GRBType.ADP, aligner=aligner)
aligner.load_reference_data(reference_loader)

# EXECUTE THE PREDICTOR
aligner.predictor()

#Get GeoJson-object of resulting prediction
featureclasses = aligner.get_results_as_geojson(resulttype=AlignerResultType.PREDICTIONS)
print (featureclasses['result'])
# SAVE PREDICTIONS TO file (GeoJSON)
aligner.save_results(path = 'output/',resulttype=AlignerResultType.PREDICTIONS)


# if fcs is None or "result" not in fcs:
#     print("empty predictions")
# else:
#     print(fcs["result"])
#     for key in dict_predictions:
#         plot_dict_diffs({key: diffs[key]})
#         show_map(
#             {key: dict_predictions[key]},
#             {key: aligner.dict_thematic[key]},
#             aligner.dict_reference,
#         )
