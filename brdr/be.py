#TODO - cleanup full file: make generic for all types, extract WFS - GML function, CRS ...
import json
from datetime import datetime
from io import BytesIO

import geopandas as gpd
import requests

from brdr.constants import (
    DATE_FORMAT,
    BRK_MAX_REFERENCE_BUFFER,
    VERSION_DATE,
    DOWNLOAD_LIMIT,
)
from brdr.geometry_utils import buffer_pos, get_bbox
from brdr.loader import GeoJsonLoader
from brdr.utils import get_collection_by_partition


def gml_response_to_geojson(url):
    # Haal de GML-response op
    response = requests.get(url)
    response.raise_for_status()  # geeft fout als de request mislukt

    # Lees de GML rechtstreeks uit de response
    gdf = gpd.read_file(BytesIO(response.content))

    # Zet om naar GeoJSON (als dict)
    return json.loads(gdf.to_json())


# https://ccff02.minfin.fgov.be/geoservices/arcgis/services/WMS/Cadastral_LayersWFS/MapServer/WFSServer?SERVICE=WFS&REQUEST=GetFeature&VERSION=2.0.0&TYPENAMES=CL:Cadastral_parcel&SRSNAME=urn:ogc:def:crs:EPSG::3812&BBOX=673571.04613103601150215,670958.87597115384414792,674256.50507965590804815,671696.05491833412088454,urn:ogc:def:crs:EPSG::3812
def get_collection_cadastral(
    geometry,
    crs="EPSG:3812"
):
    CADASTRAL_CRS = "EPSG:3812"
    name_reference_id = "CaPaKey"
    wfs = "http://ccff02.minfin.fgov.be/geoservices/arcgis/services/WMS/Cadastral_LayersWFS/MapServer/WFSServer?request=GetCapabilities&service=WFS"

    if crs == CADASTRAL_CRS:
        crs = "urn:ogc:def:crs:EPSG::3812"
    else:
        raise ValueError (f"CRS expected: {CADASTRAL_CRS}, got CRS {crs} instead")
    bbox=get_bbox(geometry)
    url = f"https://ccff02.minfin.fgov.be/geoservices/arcgis/services/WMS/Cadastral_LayersWFS/MapServer/WFSServer?SERVICE=WFS&REQUEST=GetFeature&VERSION=2.0.0&TYPENAMES=CL:Cadastral_parcel&SRSNAME={crs}&BBOX={bbox},{crs}"

    geojson = gml_response_to_geojson(url)
    collection=geojson
    return collection, name_reference_id

class BeCadastralParcelLoader(GeoJsonLoader):
    def __init__(self, aligner, partition: int = 1000):
        super().__init__()
        self.aligner = aligner
        self.part = partition
        self.data_dict_source["source"] = "Kadaster"
        #self.versiondate_info = {"name": "LastUpdDTS", "format": DATE_FORMAT}

    def load_data(self):
        if not self.aligner.dict_thematic:
            raise ValueError("Thematic data not loaded")
        # TODO CRS support for other?
        # CRS:
        # EPSG:31370 – Belge Lambert 72
        # EPSG:3812 – Belgian Lambert 2008
        # EPSG:4258 – ETRS89 (geografische coördinaten)
        supported_crs = ["EPSG:31370","EPSG:3812","EPSG:4258"]
        if self.aligner.CRS not in supported_crs:
            raise ValueError(f"BRKLoader only supports alignment in CRS '{str(supported_crs)}' while CRS '{self.aligner.CRS}' is used")
        geom_union = buffer_pos(
            self.aligner.get_thematic_union(), BRK_MAX_REFERENCE_BUFFER
        )

        collection, id_property = get_collection_cadastral(
            geometry=geom_union,
            crs=self.aligner.CRS,
        )
        self.id_property = id_property
        self.input = dict(collection)
        self.data_dict_source[VERSION_DATE] = datetime.now().strftime(DATE_FORMAT)
        self.aligner.logger.feedback_info(f"Cadaster downloaded")
        return super().load_data()
