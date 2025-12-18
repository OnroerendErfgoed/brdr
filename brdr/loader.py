"""
Core loading infrastructure for the brdr engine.

This module provides the base `Loader` abstract class and various implementations
for loading spatial data from dictionaries, GeoJSON files, URLs, and OGC services.
"""

import json
import uuid
import xml.etree.ElementTree as ET
from abc import ABC
from datetime import datetime
from typing import Any, Optional, Dict, Union

import requests
from shapely import force_2d, make_valid
from shapely.geometry.base import BaseGeometry

from brdr.constants import DATE_FORMAT, DOWNLOAD_LIMIT, MAX_REFERENCE_BUFFER, VERSION_DATE
from brdr.feature_data import AlignerFeature, AlignerFeatureCollection
from brdr.geometry_utils import buffer_pos, from_crs, to_crs
from brdr.typings import FeatureCollection, ThematicId
from brdr.utils import geojson_to_dicts, get_collection_by_partition


class Loader(ABC):
    """
    Abstract base class for all data loaders.

    The Loader handles the conversion of raw spatial data into an
    `AlignerFeatureCollection`, ensuring geometries are valid and 2D.

    Attributes:
        data_dict (Dict[ThematicId, BaseGeometry]): Mapping of IDs to geometries.
        data_dict_properties (Dict[ThematicId, dict]): Mapping of IDs to property dictionaries.
        data_dict_source (Dict[Any, str]): Metadata regarding the data source.
        versiondate_info (Optional[Dict[Any, str]]): Configuration for parsing version dates.
        is_reference (bool): Whether the loaded data serves as a reference layer.
    """

    def __init__(self, is_reference: bool):
        """
        Initialize the base Loader.

        Args:
            is_reference: Set to True if this data is a reference (target) layer
                for alignment.
        """
        self.data_dict: dict[ThematicId, BaseGeometry] = {}
        self.data_dict_properties: dict[ThematicId, dict] = {}
        self.data_dict_source: dict[Any, str] = {}
        self.versiondate_info: Optional[dict[Any, str]] = None
        self.is_reference = is_reference

    def load_data(self) -> AlignerFeatureCollection:
        """
        Processes raw data into a validated AlignerFeatureCollection.

        This method performs the following steps:
        1. Forces geometries to 2D and fixes invalid geometries.
        2. Extracts and formats version dates based on `versiondate_info`.
        3. Wraps data into `AlignerFeature` objects with unique IDs.

        Returns:
            A collection of processed features ready for alignment.
        """
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
                    except Exception:
                        # Fallback: extract only the first 10 characters (YYYY-MM-DD)
                        date = datetime.strptime(
                            self.data_dict_properties[key][
                                self.versiondate_info["name"]
                            ][:10],
                            DATE_FORMAT,
                        )

                    self.data_dict_properties[key][VERSION_DATE] = datetime.strftime(
                        date, DATE_FORMAT
                    )
                except Exception:
                    pass

        features = {
            key: AlignerFeature(
                brdr_id=uuid.uuid4().hex,
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
    """
    Loader that accepts data directly as Python dictionaries.
    """

    def __init__(
            self,
            data_dict: Dict[str, BaseGeometry],
            data_dict_properties: Dict[str, dict] = {},
            data_dict_source: Dict[str, str] = {},
            is_reference: bool = False,
    ):
        """
        Args:
            data_dict: Dictionary mapping IDs to Shapely geometries.
            data_dict_properties: Metadata properties for each feature.
            data_dict_source: Metadata about the source of the data.
            is_reference: Whether this is a reference layer.
        """
        super().__init__(is_reference=is_reference)
        self.data_dict = data_dict
        self.data_dict_properties = data_dict_properties
        self.data_dict_source = data_dict_source

    def load_data(self) -> AlignerFeatureCollection:
        """Executes the standard loading process."""
        return super().load_data()


class GeoJsonLoader(Loader):
    """
    Base class for loaders dealing with GeoJSON-formatted data.
    """

    def __init__(
            self,
            *,
            id_property: Optional[str] = None,
            _input: Optional[FeatureCollection] = None,
            is_reference: bool = False,
    ):
        """
        Args:
            id_property: The GeoJSON property to use as the unique identifier.
            _input: A dictionary representing a GeoJSON FeatureCollection.
            is_reference: Whether this is a reference layer.
        """
        super().__init__(is_reference=is_reference)
        self.id_property = id_property
        self.input = _input

    def load_data(self) -> AlignerFeatureCollection:
        """Parses the GeoJSON input and returns the feature collection."""
        self._load_geojson_data()
        return super().load_data()

    def _load_geojson_data(self) -> None:
        """
        Internal method to convert GeoJSON input into internal dictionary formats.
        """
        self.data_dict, self.data_dict_properties = geojson_to_dicts(
            self.input, self.id_property
        )


class GeoJsonFileLoader(GeoJsonLoader):
    """Loads GeoJSON data from a local file."""

    def __init__(self, path_to_file: str, id_property: str, is_reference: bool = False):
        """
        Args:
            path_to_file: System path to the .geojson file.
            id_property: The GeoJSON property to use as the unique identifier.
            is_reference: Whether this is a reference layer.
        """
        with open(path_to_file, "r") as f:
            _input = json.load(f)
        super().__init__(
            _input=_input,
            id_property=id_property,
            is_reference=is_reference,
        )


class GeoJsonUrlLoader(GeoJsonLoader):
    """Loads GeoJSON data from a remote URL."""

    def __init__(self, url: str, id_property: str, is_reference: bool = False):
        """
        Args:
            url: The URL pointing to a GeoJSON resource.
            id_property: The GeoJSON property to use as the unique identifier.
            is_reference: Whether this is a reference layer.
        """
        _input = requests.get(url).json()
        super().__init__(
            _input=_input,
            id_property=id_property,
            is_reference=is_reference,
        )


class OGCFeatureAPIReferenceLoader(GeoJsonLoader):
    """
    Loader for OGC API - Features services.

    Fetches data using spatial filters and supports partitioned downloads.
    """

    def __init__(
            self,
            url: str,
            id_property: str,
            collection: str,
            aligner: Any,
            partition: int = 1000,
            limit: int = DOWNLOAD_LIMIT,
            is_reference: bool = False,
    ):
        """
        Args:
            url: Base URL of the OGC API service.
            id_property: Property to use as unique identifier.
            collection: The specific collection ID to fetch.
            aligner: The Aligner instance for spatial context.
            partition: Number of features per request.
            limit: Maximum total features to download.
            is_reference: Whether this is a reference layer.
        """
        super().__init__(is_reference=is_reference)
        self.aligner = aligner
        self.url = url
        self.id_property = id_property
        self.part = partition
        self.limit = limit
        self.coll = collection
        self.data_dict_source["source"] = url

    def load_data(self) -> AlignerFeatureCollection:
        """
        Validates the OGC collection and CRS, then downloads data within
        the aligner's spatial extent.

        Raises:
            ValueError: If thematic data is missing, collection is not found,
                or CRS is unsupported.
        """
        if not self.aligner.thematic_data:
            raise ValueError("Thematic data not loaded")

        collections_url = self.url.rstrip("/") + "/collections"
        response = requests.get(collections_url)
        data = response.json()
        collections = [c["id"] for c in data.get("collections", [])]

        if self.coll not in collections:
            raise ValueError(f"Collection {self.coll} not found in {self.url}")

        collection_url = f"{collections_url}/{self.coll}"
        response = requests.get(collection_url)
        data = response.json()

        supported_crs = {to_crs(crs) for crs in data.get("crs", [])}
        if self.aligner.crs not in supported_crs:
            raise ValueError(f"Unsupported CRS. Supported: {supported_crs}")

        geom_union = buffer_pos(self.aligner.thematic_data.union, MAX_REFERENCE_BUFFER)
        ogcfeature_url = f"{collection_url}/items?"
        params = {"limit": self.limit, "crs": from_crs(self.aligner.crs), "f": "json"}

        collection = get_collection_by_partition(
            url=ogcfeature_url,
            params=params,
            geometry=geom_union,
            partition=self.part,
            crs=self.aligner.crs,
        )
        self.input = dict(collection)
        self.data_dict_source[VERSION_DATE] = datetime.now().strftime(DATE_FORMAT)
        self.aligner.logger.feedback_info(f"Downloaded from OGC API: {collection_url}")
        return super().load_data()


class WFSReferenceLoader(GeoJsonLoader):
    """
    Loader for Web Feature Service (WFS) version 2.0.0.
    """

    def __init__(
            self,
            url: str,
            id_property: str,
            typename: str,
            aligner: Any,
            partition: int = 1000,
            limit: int = DOWNLOAD_LIMIT,
            is_reference: bool = False,
    ):
        """
        Args:
            url: Base URL of the WFS service.
            id_property: Property to use as identifier.
            typename: The layer/feature type name.
            aligner: The Aligner instance for spatial context.
            partition: Number of features per request.
            limit: Maximum features to download.
            is_reference: Whether this is a reference layer.
        """
        super().__init__(is_reference=is_reference)
        self.aligner = aligner
        self.url = url
        self.id_property = id_property
        self.part = partition
        self.typename = typename
        self.data_dict_source["source"] = url
        self.limit = limit

    def load_data(self) -> AlignerFeatureCollection:
        """
        Parses WFS capabilities, validates CRS, and downloads features via GetFeature.

        Raises:
            ValueError: If thematic data is missing, typename is not found,
                or CRS is unsupported.
        """
        if not self.aligner.thematic_data:
            raise ValueError("Thematic data not loaded")

        params = {"service": "WFS", "version": "2.0.0", "request": "GetCapabilities"}
        response = requests.get(self.url, params=params)
        root = ET.fromstring(response.content)

        typename_exists = False
        supported_crs = []
        for feature_type in root.findall(".//{*}FeatureType"):
            name_node = feature_type.find("{*}Name")
            if name_node is not None and name_node.text == self.typename:
                typename_exists = True
                default_crs = feature_type.find("{*}DefaultCRS")
                if default_crs is not None:
                    supported_crs.append(to_crs(default_crs.text))
                for ocrs in feature_type.findall("{*}OtherCRS"):
                    try: supported_crs.append(to_crs(ocrs.text))
                    except: pass
                break

        if not typename_exists:
            raise ValueError(f"Typename {self.typename} not found in {self.url}")

        if self.aligner.crs not in supported_crs:
            raise ValueError(f"Unsupported CRS. Supported: {supported_crs}")

        geom_union = buffer_pos(self.aligner.thematic_data.union, MAX_REFERENCE_BUFFER)
        params = {
            "SERVICE": "WFS",
            "REQUEST": "GetFeature",
            "VERSION": "2.0.0",
            "TYPENAMES": self.typename,
            "SRSNAME": from_crs(self.aligner.crs, format="epsg"),
            "outputFormat": "application/json",
            "limit": self.limit,
        }

        collection = get_collection_by_partition(
            url=self.url,
            params=params,
            geometry=geom_union,
            partition=self.part,
            crs=self.aligner.crs,
        )
        self.input = dict(collection)
        self.data_dict_source[VERSION_DATE] = datetime.now().strftime(DATE_FORMAT)
        self.aligner.logger.feedback_info(f"Downloaded from WFS: {self.url}")
        return super().load_data()