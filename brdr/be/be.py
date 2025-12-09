# TODO - cleanup full file: make generic for all types, extract WFS - GML function, CRS ...
import json
from datetime import datetime
from io import BytesIO

import geopandas as gpd
import requests

from brdr.constants import (
    DATE_FORMAT,
    VERSION_DATE,
    MAX_REFERENCE_BUFFER,
)
from brdr.geometry_utils import buffer_pos, get_bbox, from_crs, to_crs
from brdr.loader import GeoJsonLoader


def gml_response_to_geojson(url, params):
    # Haal de GML-response op
    response = requests.get(url, params)
    response.raise_for_status()  # geeft fout als de request mislukt

    # Lees de GML rechtstreeks uit de response
    gdf = gpd.read_file(BytesIO(response.content))

    # Zet om naar GeoJSON (als dict)
    return json.loads(gdf.to_json())


# https://ccff02.minfin.fgov.be/geoservices/arcgis/services/WMS/Cadastral_LayersWFS/MapServer/WFSServer?SERVICE=WFS&REQUEST=GetFeature&VERSION=2.0.0&TYPENAMES=CL:Cadastral_parcel&SRSNAME=urn:ogc:def:crs:EPSG::3812&BBOX=673571.04613103601150215,670958.87597115384414792,674256.50507965590804815,671696.05491833412088454,urn:ogc:def:crs:EPSG::3812
def get_collection_cadastral(geometry, crs="EPSG:3812"):
    crs = to_crs(crs)
    name_reference_id = "CaPaKey"
    bbox = get_bbox(geometry)
    url = f"https://ccff02.minfin.fgov.be/geoservices/arcgis/services/WMS/Cadastral_LayersWFS/MapServer/WFSServer?"

    params = {
        "SERVICE": "WFS",
        "REQUEST": "GetFeature",
        "VERSION": "2.0.0",
        "TYPENAMES": "CL:Cadastral_parcel",
        "SRSNAME": from_crs(crs),
        "BBOX": bbox + "," + from_crs(crs),
        # "outputFormat": "application/json" #not available in json for this WFS, only returning XML/GML
    }
    geojson = gml_response_to_geojson(url, params)
    collection = geojson
    return collection, name_reference_id


class BeCadastralParcelLoader(GeoJsonLoader):
    def __init__(self, aligner, partition: int = 1000):
        super().__init__()
        self.aligner = aligner
        self.part = partition
        self.data_dict_source["source"] = "Kadaster"
        # self.versiondate_info = {"name": "LastUpdDTS", "format": DATE_FORMAT}

    def load_data(self):
        if not self.aligner.dict_thematic:
            raise ValueError("Thematic data not loaded")
        # CRS:
        # EPSG:31370 – Belge Lambert 72
        # EPSG:3812 – Belgian Lambert 2008
        # EPSG:4258 – ETRS89 (geografische coördinaten)
        supported_crs = ["EPSG:31370", "EPSG:3812", "EPSG:4258"]
        aligner_crs_epsg = from_crs(self.aligner.CRS, format="epsg")
        if aligner_crs_epsg not in supported_crs:
            raise ValueError(
                f"BeCadastralParcelLoader only supports alignment in CRS '{str(supported_crs)}' while CRS '{aligner_crs_epsg}' is used"
            )
        geom_union = buffer_pos(self.aligner.get_thematic_union(), MAX_REFERENCE_BUFFER)

        collection, id_property = get_collection_cadastral(
            geometry=geom_union,
            crs=self.aligner.CRS,
        )
        self.id_property = id_property
        self.input = dict(collection)
        self.data_dict_source[VERSION_DATE] = datetime.now().strftime(DATE_FORMAT)
        self.aligner.logger.feedback_info(f"Cadaster downloaded")
        return super().load_data()
