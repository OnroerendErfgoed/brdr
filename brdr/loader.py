"""
Core loading infrastructure for the brdr engine.

This module provides the base `Loader` abstract class and various implementations
for loading spatial data from dictionaries, GeoJSON files, URLs, and OGC services.
"""

import json
import os
import uuid
import xml.etree.ElementTree as ET
from abc import ABC
from datetime import datetime
from typing import Any, Optional, Dict

import requests
from shapely import force_2d, make_valid
from shapely.geometry.base import BaseGeometry

from brdr.constants import (
    DATE_FORMAT,
    DOWNLOAD_LIMIT,
    MAX_REFERENCE_BUFFER,
    VERSION_DATE,
)
from brdr.feature_data import AlignerFeature, AlignerFeatureCollection
from brdr.geometry_utils import buffer_pos, from_crs, to_crs
from brdr.typings import FeatureCollection, InputId
from brdr.utils import geojson_to_dicts, get_collection_by_partition


class Loader(ABC):
    """
    Abstract base class for all data loaders.

    The Loader handles the conversion of raw spatial data into an
    `AlignerFeatureCollection`, ensuring geometries are valid and 2D.

    Attributes
    ----------
    data_dict : Dict[ThematicId, BaseGeometry]
        Mapping of unique IDs to their corresponding Shapely geometries.
    data_dict_properties : Dict[ThematicId, dict]
        Mapping of IDs to attribute dictionaries.
    data_dict_source : Dict[Any, str]
        Metadata regarding the data source (e.g., URL or file path).
    versiondate_info : dict, optional
        Configuration for parsing version dates, containing 'name' (field)
        and 'format' (strptime pattern).
    is_reference : bool
        Whether the loaded data serves as a reference (target) layer.

    Notes
    -----
    All subclasses must implement `load_data` or populate the internal
    data dictionaries before calling the base `load_data` method.
    """

    def __init__(self, is_reference: bool):
        """
        Initializes the base Loader.

        Parameters
        ----------
        is_reference : bool
            Set to True if this data is a reference layer for alignment.
        """
        self.data_dict: dict[InputId, BaseGeometry] = {}
        self.data_dict_properties: dict[InputId, dict] = {}
        self.data_dict_source: dict[Any, str] = {
            "source_url": None,
        }
        self.versiondate_info: Optional[dict[Any, str]] = None
        self.is_reference = is_reference

    def load_data(self) -> AlignerFeatureCollection:
        """
        Processes raw data into a validated AlignerFeatureCollection.

        This method performs geometric validation, version date extraction,
        and wraps the data into AlignerFeature objects.

        Returns
        -------
        AlignerFeatureCollection
            A collection of processed and validated features.

        Notes
        -----
        The processing pipeline follows these steps:



        1. **Geometry Fix**: Forces 2D and applies `make_valid`.
        2. **Metadata Enrichment**: Parses version dates if `versiondate_info` is set.
        3. **Feature Creation**: Generates unique UUIDs and wraps features.
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
                data_id=key,
                brdr_id=uuid.uuid4().urn,
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
    Loader that accepts data directly from Python dictionaries.

    Useful for integrating with existing in-memory data structures.
    """

    def __init__(
        self,
        data_dict: Dict[str, BaseGeometry],
        data_dict_properties: Dict[str, dict] = {},
        data_dict_source: Dict[str, str] = {"source_url": None},
        is_reference: bool = False,
    ):
        """
        Parameters
        ----------
        data_dict : Dict[str, BaseGeometry]
            Dictionary mapping IDs to Shapely geometries.
        data_dict_properties : Dict[str, dict], optional
            Metadata properties for each feature.
        data_dict_source : Dict[str, str], optional
            Metadata about the data origin.
        is_reference : bool, optional
            Whether this is a reference layer. Defaults to False.
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

    Processes a GeoJSON FeatureCollection into the internal dictionary format.
    """

    def __init__(
        self,
        *,
        id_property: Optional[str] = None,
        _input: Optional[FeatureCollection] = None,
        source: Optional[str] = None,
        source_url: Optional[str] = None,
        is_reference: bool = False,
    ):
        """
        Parameters
        ----------
        id_property : str, optional
            The GeoJSON property key to use as the unique identifier.
        _input : FeatureCollection, optional
            A dictionary following the GeoJSON FeatureCollection spec.
        is_reference : bool, optional
            Whether this is a reference layer. Defaults to False.
        """
        super().__init__(is_reference=is_reference)
        self.id_property = id_property
        self.input = _input
        self.data_dict_source["source"] = source
        self.data_dict_source["source_url"] = source_url

    def load_data(self) -> AlignerFeatureCollection:
        """Parses the GeoJSON input and returns the feature collection."""
        self._load_geojson_data()
        return super().load_data()

    def _load_geojson_data(self) -> None:
        """Internal method to convert GeoJSON input into internal dictionary formats."""
        self.data_dict, self.data_dict_properties = geojson_to_dicts(
            self.input, self.id_property
        )


class GeoJsonFileLoader(GeoJsonLoader):
    """Loads GeoJSON data from a local file system."""

    def __init__(self, path_to_file: str, id_property: str, is_reference: bool = False):
        """
        Parameters
        ----------
        path_to_file : str
            Full system path to the `.geojson` file.
        id_property : str
            The GeoJSON property key to use as the identifier.
        is_reference : bool, optional
            Whether this is a reference layer. Defaults to False.
        """
        with open(path_to_file, "r") as f:
            _input = json.load(f)
        super().__init__(
            _input=_input,
            source=path_to_file,
            source_url=f"file://{os.path.abspath(path_to_file)}",
            id_property=id_property,
            is_reference=is_reference,
        )


class GeoJsonUrlLoader(GeoJsonLoader):
    """Loads GeoJSON data from a remote URL via HTTP GET."""

    def __init__(self, url: str, id_property: str, is_reference: bool = False):
        """
        Parameters
        ----------
        url : str
            The URL pointing to a GeoJSON resource.
        id_property : str
            The GeoJSON property key to use as the identifier.
        is_reference : bool, optional
            Whether this is a reference layer. Defaults to False.
        """
        _input = requests.get(url).json()
        super().__init__(
            _input=_input,
            source=url,
            source_url=url,
            id_property=id_property,
            is_reference=is_reference,
        )


class OGCFeatureAPIReferenceLoader(GeoJsonLoader):
    """
    Loader for OGC API - Features services.

    Fetches data using spatial filters (BBOX) and handles partitioning
    of large datasets based on the thematic extent.
    """

    def __init__(
        self,
        url: str,
        id_property: str,
        collection: str,
        aligner: Any,
        partition: int = 1000,
        limit: int = DOWNLOAD_LIMIT,
        is_reference: bool = True,
    ):
        """
        Parameters
        ----------
        url : str
            Base URL of the OGC API service.
        id_property : str
            Property to use as unique identifier.
        collection : str
            The specific collection ID to fetch from the service.
        aligner : Aligner
            The Aligner instance, used for CRS and spatial extent context.
        partition : int, optional
            Number of features to request per network call. Defaults to 1000.
        limit : int, optional
            Maximum total features to download. Defaults to DOWNLOAD_LIMIT.
        is_reference : bool, optional
            Whether this is a reference layer: True.

        """
        super().__init__(is_reference=is_reference)
        self.aligner = aligner
        self.url = url
        self.id_property = id_property
        self.part = partition
        self.limit = limit
        self.coll = collection
        self.data_dict_source["source"] = url
        self.data_dict_source["source_url"] = url + "/" + collection

    def load_data(self) -> AlignerFeatureCollection:
        """
        Validates the OGC service and downloads data within the thematic extent.

        Returns
        -------
        AlignerFeatureCollection
            The downloaded and processed reference data.

        Raises
        ------
        ValueError
            If thematic data is missing, the collection is not found,
            or the service does not support the Aligner's CRS.

        Notes
        -----

        ```mermaid
        graph LR
        Thematic[Thematic Data Extent] --> Buffer[Apply MAX_REFERENCE_BUFFER]
        Buffer --> Partition[Split into BBOX Partitions]
        Partition --> Request[HTTP GetFeature/Items]
        Request --> Merge[Merge Results into Collection]

        subgraph Service[External OGC/WFS Service]
        Request
        end
        ```

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

        supported_crs = []
        for crs in data.get("crs", "0000"):
            try:
                supported_crs.append(to_crs(crs))
            except:
                pass
        if self.aligner.crs not in supported_crs:
            raise ValueError(f"Unsupported CRS. Supported: {supported_crs}")

        geom_union = buffer_pos(self.aligner.thematic_data.union, MAX_REFERENCE_BUFFER)
        if geom_union is None or geom_union.is_empty:
            raise ValueError(
                "Reference could not be loaded. Please load thematic data first"
            )
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
    Loader for OGC Web Feature Service (WFS) version 2.0.0.

    Supports GetCapabilities parsing and GetFeature requests with spatial
    filtering using the Aligner's extent.
    """

    def __init__(
        self,
        url: str,
        id_property: str,
        typename: str,
        aligner: Any,
        partition: int = 1000,
        limit: int = DOWNLOAD_LIMIT,
        is_reference: bool = True,
    ):
        """
        Parameters
        ----------
        url : str
            Base URL of the WFS service.
        id_property : str
            Property to use as unique identifier.
        typename : str
            The layer/feature type name (e.g., 'ns:layer_name').
        aligner : Aligner
            The Aligner instance for spatial context.
        partition : int, optional
            Number of features per request. Defaults to 1000.
        limit : int, optional
            Maximum total features to download. Defaults to DOWNLOAD_LIMIT.
        is_reference : bool, optional
            Whether this is a reference layer: True.
        """
        super().__init__(is_reference=is_reference)
        self.aligner = aligner
        self.url = url
        self.id_property = id_property
        self.part = partition
        self.typename = typename
        self.data_dict_source["source"] = url
        self.data_dict_source["source_url"] = url
        self.limit = limit

    def load_data(self) -> AlignerFeatureCollection:
        """
        Downloads features via GetFeature based on capabilities and extent.

        Returns
        -------
        AlignerFeatureCollection
            The downloaded and processed reference data.

        Raises
        ------
        ValueError
            If thematic data is missing, typename is not found,
            or the service does not support the Aligner's CRS.
        """
        if not self.aligner.thematic_data:
            raise ValueError("Thematic data not loaded")
        self.data_dict_source[VERSION_DATE] = datetime.now().strftime(DATE_FORMAT)
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
                    try:
                        supported_crs.append(to_crs(ocrs.text))
                    except:
                        pass
                break

        if not typename_exists:
            raise ValueError(f"Typename {self.typename} not found in {self.url}")

        if self.aligner.crs not in supported_crs:
            raise ValueError(f"Unsupported CRS. Supported: {supported_crs}")

        geom_union = buffer_pos(self.aligner.thematic_data.union, MAX_REFERENCE_BUFFER)
        if geom_union is None or geom_union.is_empty:
            raise ValueError(
                "Reference could not be loaded. Please load thematic data first"
            )
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
        self.aligner.logger.feedback_info(f"Downloaded from OGC WFS: {self.url}")
        return super().load_data()
