import json
import uuid
import xml.etree.ElementTree as ET
from abc import ABC
from datetime import datetime
from typing import Any

import requests as requests
from shapely import force_2d
from shapely import make_valid
from shapely.geometry.base import BaseGeometry

from brdr.constants import DATE_FORMAT
from brdr.constants import DOWNLOAD_LIMIT
from brdr.constants import MAX_REFERENCE_BUFFER
from brdr.constants import VERSION_DATE
from brdr.feature_data import AlignerFeature
from brdr.feature_data import AlignerFeatureCollection
from brdr.geometry_utils import buffer_pos
from brdr.geometry_utils import from_crs
from brdr.geometry_utils import to_crs
from brdr.typings import FeatureCollection
from brdr.typings import ThematicId
from brdr.utils import geojson_to_dicts
from brdr.utils import get_collection_by_partition


class Loader(ABC):
    def __init__(self, is_reference: bool):
        self.data_dict: dict[ThematicId, BaseGeometry] = {}
        self.data_dict_properties: dict[ThematicId, dict] = {}
        self.data_dict_source: dict[Any, str] = {}
        self.versiondate_info: dict[Any, str] = None
        self.is_reference = is_reference

    def load_data(self):
        self.data_dict = {
            x: force_2d(make_valid(self.data_dict[x])) for x in self.data_dict
        }
        if self.versiondate_info is not None:
            for key in self.data_dict_properties.keys():
                try:
                    try:
                        date = datetime.strptime(
                            self.data_dict_properties[key][
                                self.versiondate_info["name"]
                            ],
                            self.versiondate_info["format"],
                        )
                    except:
                        # Catch, to try extracting only the date with default -date format if specific format does not work

                        date = datetime.strptime(
                            self.data_dict_properties[key][
                                self.versiondate_info["name"]
                            ][:10],
                            DATE_FORMAT,
                        )

                    self.data_dict_properties[key][VERSION_DATE] = datetime.strftime(
                        date, DATE_FORMAT
                    )
                except:
                    # No version date added to features
                    pass

        return self.data_dict, self.data_dict_properties, self.data_dict_source

    def load_data_as_feature_collection(self):
        # TODO: rework loaders to make load_data itself return AlignerFeatureCollection
        self.load_data()
        features = {
            key: AlignerFeature(
                brdr_id=uuid.uuid4().hex, #TODO use id from dict if available
                geometry=self.data_dict[key],
                properties=self.data_dict_properties.get(key, {}),
            )
            for key in self.data_dict
        }
        return AlignerFeatureCollection(
            features=features,
            source=self.data_dict_source,
            is_reference=self.is_reference,
        )


class DictLoader(Loader):
    def __init__(
        self,
        data_dict: dict[str:BaseGeometry],
        data_dict_properties: dict[str:dict] = {},
        data_dict_source: dict[str:str] = {},
        is_reference: bool = False,
    ):
        super().__init__(is_reference=is_reference)
        self.data_dict = data_dict
        self.data_dict_properties = data_dict_properties
        self.data_dict_source = data_dict_source

    def load_data(self):
        return super().load_data()


class GeoJsonLoader(Loader):
    def __init__(
        self,
        *,
        id_property: str = None,
        _input: FeatureCollection = None,
        is_reference: bool = False,
    ):
        super().__init__(is_reference=is_reference)
        self.id_property = id_property
        self.input = _input

    def load_data(self):
        self._load_geojson_data()
        return super().load_data()

    def _load_geojson_data(self):
        """
        Load geometries of a GeoJSON and stores them in a dictionary.

        This method processes the thematic data from the input GeoJSON file. It
        iterates through each feature, extracts the relevant properties, converts the
        geometry to a valid shape, and stores it in a dictionary.

        Returns:
            None.
        """
        # THEMATIC PREPARATION
        self.data_dict, self.data_dict_properties = geojson_to_dicts(
            self.input, self.id_property
        )
        return


class GeoJsonFileLoader(GeoJsonLoader):
    def __init__(
        self,
        path_to_file,
        id_property,
        is_reference: bool = False,
    ):
        with open(path_to_file, "r") as f:
            _input = json.load(f)
        super().__init__(
            _input=_input,
            id_property=id_property,
            is_reference=is_reference,
        )


class GeoJsonUrlLoader(GeoJsonLoader):

    def __init__(
        self,
        url,
        id_property,
        is_reference: bool = False,
    ):
        _input = requests.get(url).json()
        super().__init__(
            _input=_input,
            id_property=id_property,
            is_reference=is_reference,
        )


class OGCFeatureAPIReferenceLoader(GeoJsonLoader):

    def __init__(
        self,
        url,
        id_property,
        collection,
        aligner,
        partition: int = 1000,
        limit=DOWNLOAD_LIMIT,
        is_reference: bool = False,
    ):
        super().__init__(is_reference=is_reference)
        self.aligner = aligner
        self.url = url
        self.id_property = id_property
        self.part = partition
        self.limit = limit
        self.coll = collection
        self.data_dict_source["source"] = url

    def load_data(self):
        if not self.aligner.dict_thematic:
            raise ValueError("Thematic data not loaded")
        # Check if collection exists
        collections_url = self.url.rstrip("/") + "/" + "collections"
        response = requests.get(collections_url)
        data = response.json()
        collections = [collection["id"] for collection in data.get("collections", [])]
        if self.coll not in collections:
            raise ValueError(
                f"Collection {self.coll} not found inside OGC Feature API {self.url}"
            )
        collection_url = collections_url + "/" + self.coll

        # Get all supported CRS for the collection of OGCFeatureAPI
        response = requests.get(collection_url)
        data = response.json()
        supported_crs = set()
        for crs in data.get("crs", []):
            try:
                supported_crs.add(to_crs(crs))
            except:
                pass
        if self.aligner.CRS not in supported_crs:
            raise ValueError(
                f"OGCFeatureAPIReferenceLoader '{collection_url}' only supports alignment in CRS '{str(supported_crs)}' while CRS '{self.aligner.CRS}' is used."
            )
        geom_union = buffer_pos(self.aligner.get_thematic_union(), MAX_REFERENCE_BUFFER)
        ogcfeature_url = collection_url + "/items?"
        params = {"limit": self.limit, "crs": from_crs(self.aligner.CRS), "f": "json"}

        collection = get_collection_by_partition(
            url=ogcfeature_url,
            params=params,
            geometry=geom_union,
            partition=self.part,
            crs=self.aligner.CRS,
        )
        self.input = dict(collection)
        self.data_dict_source[VERSION_DATE] = datetime.now().strftime(DATE_FORMAT)
        self.aligner.logger.feedback_info(
            f"Reference data downloaded from OGC Feature API: {collection_url}"
        )
        return super().load_data()


class WFSReferenceLoader(GeoJsonLoader):

    def __init__(
        self,
        url,
        id_property,
        typename,
        aligner,
        partition: int = 1000,
        limit=DOWNLOAD_LIMIT,
        is_reference: bool = False,
    ):
        super().__init__(is_reference=is_reference)
        self.aligner = aligner
        self.url = url
        self.id_property = id_property
        self.part = partition
        self.typename = typename
        self.data_dict_source["source"] = url
        self.limit = limit

    def load_data(self):
        if not self.aligner.dict_thematic:
            raise ValueError("Thematic data not loaded")

        # Get all supported CRS for the collection of OGCFeatureAPI

        params = {"service": "WFS", "version": "2.0.0", "request": "GetCapabilities"}

        response = requests.get(self.url, params=params)

        # Ophalen van de XML
        root = ET.fromstring(response.content)
        typename_exists = False
        for feature_type in root.findall(".//{*}FeatureType"):
            for child in feature_type:
                if child.tag.endswith("Name") and child.text == self.typename:
                    typename_exists = True
                    supported_crs = []
                    defaultcrs = feature_type.find("{*}DefaultCRS")
                    supported_crs.append(to_crs(defaultcrs.text))
                    for ocrs in feature_type.findall("{*}OtherCRS"):
                        try:
                            supported_crs.append(to_crs(ocrs.text))
                        except:
                            pass
                    break
        if not typename_exists:
            raise ValueError(
                f"Collection {self.typename} not found inside OGC WFS{self.url}"
            )

        if self.aligner.CRS not in supported_crs:
            raise ValueError(
                f"WFS '{self.url}' only supports alignment in CRS '{str(supported_crs)}' while CRS '{self.aligner.CRS}' is used."
            )
        geom_union = buffer_pos(self.aligner.get_thematic_union(), MAX_REFERENCE_BUFFER)

        params = {
            "SERVICE": "WFS",
            "REQUEST": "GetFeature",
            "VERSION": "2.0.0",
            "TYPENAMES": self.typename,
            "SRSNAME": from_crs(self.aligner.CRS, format="epsg"),
            "outputFormat": "application/json",
            "limit": self.limit,
        }

        collection = get_collection_by_partition(
            url=self.url,
            params=params,
            geometry=geom_union,
            partition=self.part,
            crs=self.aligner.CRS,
        )
        self.input = dict(collection)
        self.data_dict_source[VERSION_DATE] = datetime.now().strftime(DATE_FORMAT)
        self.aligner.logger.feedback_info(
            f"Reference data downloaded from WFS: {self.url}"
        )
        return super().load_data()
