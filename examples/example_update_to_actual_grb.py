from brdr.aligner import Aligner
from brdr.grb import GRBFiscalParcelLoader
from brdr.grb import update_to_actual_grb
from brdr.loader import GeoJsonFileLoader

#Create a featurecollection (aligned on 2022), to use for the 'update_to_actual_grb'
base_year = "2022"
base_aligner = Aligner()
name_thematic_id = "theme_identifier"
loader = GeoJsonFileLoader("themelayer.geojson", name_thematic_id)
base_aligner.load_thematic_data(loader)
base_aligner.load_reference_data(
    GRBFiscalParcelLoader(year=base_year, aligner=base_aligner)
)
base_process_result = base_aligner.process_dict_thematic(relevant_distance=2)
fcs = base_aligner.get_results_as_geojson(formula=True)
featurecollection_base_result= fcs["result"]
print (featurecollection_base_result)
#Update Featurecollection to actual version
featurecollection = update_to_actual_grb(featurecollection_base_result,base_aligner.name_thematic_id)
#Print results
for feature in featurecollection["result"]["features"]:
    print(
        feature["properties"][name_thematic_id]
        + ": "
        + feature["properties"]["evaluation"]
    )
geojson = featurecollection["result"]
print(geojson)
