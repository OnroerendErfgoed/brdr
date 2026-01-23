import json
from datetime import datetime
from io import BytesIO

import geopandas as gpd
import requests

from brdr.be.constants import BE_SUPPORTED_CRS
from brdr.constants import (
    DATE_FORMAT,
    VERSION_DATE,
    MAX_REFERENCE_BUFFER,
)
from brdr.geometry_utils import buffer_pos, get_bbox, from_crs, to_crs
from brdr.loader import GeoJsonLoader


def gml_response_to_geojson(url, params):
    """
    Fetches a GML (Geography Markup Language) response from a URL,
    reads it using GeoPandas, and converts the result to
    GeoJSON format (as a Python dictionary).

    Args:
        url (str): The URL to fetch the GML data from (e.g., a WFS endpoint).
        params (dict): The parameters to be sent with the GET request
                       (e.g., service, version, request, typeName).

    Returns:
        dict: The geographical data in GeoJSON format.

    Raises:
        requests.exceptions.HTTPError: If the HTTP request fails.
        fiona.errors.DriverError: If the GML data cannot be read correctly.
    """
    # Fetch the GML response
    response = requests.get(url, params)
    response.raise_for_status()  # Raises an error if the request failed

    # Read the GML directly from the response (using BytesIO for in-memory processing)
    gdf = gpd.read_file(BytesIO(response.content))

    # Convert to GeoJSON (as a dict)
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
        if not self.aligner.crs in (to_crs(element) for element in BE_SUPPORTED_CRS):
            raise ValueError(
                f"BeCadastralParcelLoader only supports alignment in CRS '{BE_SUPPORTED_CRS}' while CRS '{self.aligner.crs}' is used"
            )
        self.part = partition
        self.data_dict_source["source"] = "Kadaster"
        # self.versiondate_info = {"name": "LastUpdDTS", "format": DATE_FORMAT}

    def load_data(self):
        if not self.aligner.thematic_data:
            raise ValueError("Thematic data not loaded")
        geom_union = buffer_pos(self.aligner.thematic_data.union, MAX_REFERENCE_BUFFER)
        if geom_union is None or geom_union.is_empty:
            raise ValueError(
                "Reference could not be loaded. Please load thematic data first"
            )
        collection, id_property = get_collection_cadastral(
            geometry=geom_union,
            crs=self.aligner.crs,
        )
        self.id_property = id_property
        self.input = dict(collection)
        self.data_dict_source[VERSION_DATE] = datetime.now().strftime(DATE_FORMAT)
        self.aligner.logger.feedback_info(f"Cadaster downloaded")
        return super().load_data()
