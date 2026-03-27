from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from typing import Optional

import requests

from brdr.be.constants import BE_SUPPORTED_CRS
from brdr.constants import (
    DATE_FORMAT,
    VERSION_DATE,
    MAX_REFERENCE_BUFFER,
    DOWNLOAD_LIMIT,
)
from brdr.geometry_utils import (
    buffer_pos,
    from_crs,
    get_bbox,
    get_partitions,
    to_crs,
)
from brdr.loader import GeoJsonLoader
from brdr.utils import deduplicate_features, get_collection_by_partition


# https://ccff02.minfin.fgov.be/geoservices/arcgis/services/WMS/Cadastral_LayersWFS/MapServer/WFSServer?SERVICE=WFS&REQUEST=GetFeature&VERSION=2.0.0&TYPENAMES=CL:Cadastral_parcel&SRSNAME=urn:ogc:def:crs:EPSG::3812&BBOX=673571.04613103601150215,670958.87597115384414792,674256.50507965590804815,671696.05491833412088454,urn:ogc:def:crs:EPSG::3812
def get_collection_cadastral(geometry, crs="EPSG:3812", partition=1000):
    crs = to_crs(crs)
    name_reference_id = "CaPaKey"
    url = f"https://ccff02.minfin.fgov.be/geoservices/arcgis/services/WMS/Cadastral_LayersWFS/MapServer/WFSServer?"

    params = {
        "SERVICE": "WFS",
        "REQUEST": "GetFeature",
        "VERSION": "2.0.0",
        "TYPENAMES": "CL:Cadastral_parcel",
        "SRSNAME": from_crs(crs),
        "limit": 10000,
        # "outputFormat": "application/json" #not available in json for this WFS, only returning XML/GML
    }
    collection = get_collection_by_partition(
        url=url,
        params=params,
        geometry=geometry,
        partition=partition,
        crs=crs,
    )
    return collection, name_reference_id


def _safe_get_arcgis_max_record_count(url: str, fallback: int = 2000) -> int:
    """
    Read maxRecordCount from an ArcGIS layer metadata endpoint.
    Falls back to a conservative value when unavailable.
    """
    try:
        response = requests.get(url, params={"f": "pjson"}, timeout=30)
        response.raise_for_status()
        data = response.json()
        max_record_count = int(data.get("maxRecordCount", fallback))
        if max_record_count <= 0:
            return fallback
        return max_record_count
    except Exception:
        return fallback


def _fetch_arcgis_partition_features(
    url: str, bbox: str, wkid: int, page_size: int, request_timeout: int
):
    """
    Fetch all features for one bbox partition from ArcGIS REST /query endpoint.
    """
    features = []
    offset = 0
    previous_page_signature = None
    while True:
        params = {
            "where": "1=1",
            "geometry": bbox,
            "geometryType": "esriGeometryEnvelope",
            "spatialRel": "esriSpatialRelIntersects",
            "inSR": wkid,
            "outSR": wkid,
            "outFields": "*",
            "returnGeometry": "true",
            "resultRecordCount": page_size,
            "resultOffset": offset,
            "f": "geojson",
        }
        response = requests.get(url, params=params, timeout=request_timeout)
        response.raise_for_status()
        data = response.json()
        page_features = data.get("features", [])
        if not page_features:
            break

        features.extend(page_features)

        # Some services ignore offset and keep returning the same page; guard against loops.
        page_signature = tuple(
            f.get("id", f.get("properties", {}).get("OBJECTID")) for f in page_features
        )
        if page_signature == previous_page_signature:
            break
        previous_page_signature = page_signature

        exceeded = bool(data.get("exceededTransferLimit", False))
        if not exceeded and len(page_features) < page_size:
            break
        offset += len(page_features)

    return features


def get_collection_cadastral_arcgis(
    geometry,
    crs="EPSG:3812",
    partition=1000,
    max_workers: Optional[int] = None,
    request_timeout: int = 60,
):
    """
    Retrieve cadastral parcels using ArcGIS REST /query with automatic page-size detection.
    """
    crs = to_crs(crs)
    name_reference_id = "CaPaKey"
    wkid = crs.to_epsg() or 3812
    url = "https://ccff02.minfin.fgov.be/geoservices/arcgis/rest/services/WMS/Cadastral_Layers/MapServer/11/query"
    page_size = _safe_get_arcgis_max_record_count(url)

    if geometry is None or geometry.is_empty:
        return {"type": "FeatureCollection", "features": []}, name_reference_id

    if partition < 1:
        geoms = [geometry]
    else:
        geoms = get_partitions(geometry, partition)

    all_features = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(
                _fetch_arcgis_partition_features,
                url,
                get_bbox(g),
                wkid,
                page_size,
                request_timeout,
            ): g
            for g in geoms
        }
        for future in as_completed(futures):
            all_features.extend(future.result())

    return {
        "type": "FeatureCollection",
        "features": deduplicate_features(all_features),
    }, name_reference_id


class BeCadastralParcelLoader(GeoJsonLoader):
    """
    Reference loader for Belgian cadastral parcels.

    This loader retrieves parcel geometries that intersect the buffered thematic
    extent and exposes them as a standard `AlignerFeatureCollection`.

    Two backends are supported:
    - `arcgis` (default): ArcGIS REST query endpoint with automatic server page-size detection.
    - `wfs`: WFS 2.0.0 endpoint using the generic OGC partition/pagination flow.

    Notes
    -----
    - Only Belgian CRS definitions listed in `BE_SUPPORTED_CRS` are supported.
    - Spatial partitioning is applied to keep requests bounded and improve throughput.
    """

    def __init__(
        self,
        aligner,
        partition: int = 1000,
        service: str = "arcgis",
        max_workers: Optional[int] = None,
        request_timeout: int = 60,
    ):
        """
        Initialize a Belgian cadastral parcel loader.

        Parameters
        ----------
        aligner : Aligner
            Aligner instance containing thematic data and CRS context.
        partition : int, optional
            Grid size parameter used to partition the thematic union extent into
            smaller request windows. Use values < 1 to disable partitioning.
        service : {"arcgis", "wfs"}, optional
            Backend service type used to retrieve parcels.
            Defaults to `"arcgis"`.
        max_workers : int, optional
            Maximum number of worker threads used for partition requests in ArcGIS mode.
            Defaults to `None` (executor default).
        request_timeout : int, optional
            HTTP timeout in seconds for ArcGIS requests.
            Defaults to 60.

        Raises
        ------
        ValueError
            If the aligner CRS is not supported or `service` is invalid.
        """
        super().__init__()
        self.aligner = aligner
        if not self.aligner.crs in (to_crs(element) for element in BE_SUPPORTED_CRS):
            raise ValueError(
                f"BeCadastralParcelLoader only supports alignment in CRS '{BE_SUPPORTED_CRS}' while CRS '{self.aligner.crs}' is used"
            )
        self.part = partition
        self.service = service.lower()
        if self.service not in {"arcgis", "wfs"}:
            raise ValueError("service must be 'arcgis' or 'wfs'")
        self.max_workers = max_workers
        self.request_timeout = request_timeout
        self.data_dict_source["source"] = "Kadaster"
        # self.versiondate_info = {"name": "LastUpdDTS", "format": DATE_FORMAT}

    def load_data(self):
        """
        Load cadastral reference data for the current thematic extent.

        The thematic union is buffered with `MAX_REFERENCE_BUFFER`, then parcel
        features are fetched from the configured backend and converted to the
        internal feature collection format.

        Returns
        -------
        AlignerFeatureCollection
            Validated reference features ready for alignment.

        Raises
        ------
        ValueError
            If thematic data is missing or no valid thematic extent is available.
        """
        if not self.aligner.thematic_data:
            raise ValueError("Thematic data not loaded")
        geom_union = buffer_pos(self.aligner.thematic_data.union, MAX_REFERENCE_BUFFER)
        if geom_union is None or geom_union.is_empty:
            raise ValueError(
                "Reference could not be loaded. Please load thematic data first"
            )
        if self.service == "arcgis":
            collection, id_property = get_collection_cadastral_arcgis(
                geometry=geom_union,
                crs=self.aligner.crs,
                partition=self.part,
                max_workers=self.max_workers,
                request_timeout=self.request_timeout,
            )
        else:
            collection, id_property = get_collection_cadastral(
                geometry=geom_union,
                crs=self.aligner.crs,
                partition=self.part,
            )
        self.id_property = id_property
        self.input = dict(collection)
        self.data_dict_source[VERSION_DATE] = datetime.now().strftime(DATE_FORMAT)
        self.aligner.logger.feedback_info(
            f"Cadaster downloaded (service={self.service})"
        )
        return super().load_data()


def get_collection_vrbg(
    geometry,
    collection: str,
    crs="EPSG:3812",
    partition=1000,
    limit: int = DOWNLOAD_LIMIT,
    max_workers: Optional[int] = None,
    max_pages=float("inf"),
    request_timeout: int = 60,
):
    """
    Retrieve administrative boundaries from the VRBG OGC API Features service.
    """
    crs = to_crs(crs)
    url = (
        "https://geo.api.vlaanderen.be/VRBG/ogc/features/v1/"
        f"collections/{collection}/items"
    )
    params = {
        "f": "json",
        "crs": from_crs(crs),
        "limit": limit,
    }
    return get_collection_by_partition(
        url=url,
        params=params,
        geometry=geometry,
        partition=partition,
        crs=crs,
        max_workers=max_workers,
        max_pages=max_pages,
        request_timeout=request_timeout,
    )


class BeAdministrativeBoundaryLoader(GeoJsonLoader):
    """
    Reference loader for Belgian administrative boundaries (VRBG).

    Data is fetched from the Vlaanderen OGC API Features endpoint and clipped by
    the buffered thematic extent through BBOX partition requests.
    """

    def __init__(
        self,
        aligner,
        collection: str,
        id_property: str = "OIDN",
        partition: int = 1000,
        limit: int = DOWNLOAD_LIMIT,
        max_workers: Optional[int] = None,
        max_pages=float("inf"),
        request_timeout: int = 60,
    ):
        """
        Initialize a VRBG administrative boundary loader.

        Parameters
        ----------
        aligner : Aligner
            Aligner instance containing thematic data and CRS context.
        collection : str
            VRBG collection identifier from `/collections` (for example a boundary layer).
        id_property : str, optional
            Feature property used as unique identifier. Defaults to `"OIDN"`.
        partition : int, optional
            Grid size parameter used to partition the thematic extent.
        limit : int, optional
            Requested number of features per request page.
        max_workers : int, optional
            Maximum number of concurrent workers for partition requests.
        max_pages : int or float, optional
            Maximum number of pages per partition. Defaults to unlimited.
        request_timeout : int, optional
            HTTP request timeout in seconds. Defaults to 60.

        Raises
        ------
        ValueError
            If the aligner CRS is not supported.
        """
        super().__init__()
        self.aligner = aligner
        if not self.aligner.crs in (to_crs(element) for element in BE_SUPPORTED_CRS):
            raise ValueError(
                f"BeAdministrativeBoundaryLoader only supports alignment in CRS '{BE_SUPPORTED_CRS}' while CRS '{self.aligner.crs}' is used"
            )
        self.collection = collection
        self.id_property = id_property
        self.partition = partition
        self.limit = limit
        self.max_workers = max_workers
        self.max_pages = max_pages
        self.request_timeout = request_timeout
        self.data_dict_source["source"] = "VRBG"
        self.data_dict_source["source_url"] = (
            "https://geo.api.vlaanderen.be/VRBG/ogc/features/v1"
        )

    def load_data(self):
        """
        Load administrative boundary reference data for the thematic extent.

        Returns
        -------
        AlignerFeatureCollection
            Validated reference features ready for alignment.

        Raises
        ------
        ValueError
            If thematic data is missing or no valid thematic extent is available.
        """
        if not self.aligner.thematic_data:
            raise ValueError("Thematic data not loaded")
        geom_union = buffer_pos(self.aligner.thematic_data.union, MAX_REFERENCE_BUFFER)
        if geom_union is None or geom_union.is_empty:
            raise ValueError(
                "Reference could not be loaded. Please load thematic data first"
            )

        collection = get_collection_vrbg(
            geometry=geom_union,
            collection=self.collection,
            crs=self.aligner.crs,
            partition=self.partition,
            limit=self.limit,
            max_workers=self.max_workers,
            max_pages=self.max_pages,
            request_timeout=self.request_timeout,
        )
        self.input = dict(collection)
        self.data_dict_source[VERSION_DATE] = datetime.now().strftime(DATE_FORMAT)
        self.aligner.logger.feedback_info(
            f"VRBG boundaries downloaded (collection={self.collection})"
        )
        return super().load_data()
