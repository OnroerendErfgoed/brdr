from datetime import datetime
from typing import Any
from typing import Iterable

import geopandas as gpd
import osmnx as ox

from brdr.constants import (
    DATE_FORMAT,
    VERSION_DATE,
)
from brdr.geometry_utils import buffer_pos
from brdr.loader import DictLoader
from brdr.osm.constants import OSM_MAX_REFERENCE_BUFFER


class OSMLoader(DictLoader):
    """
    Loader for OpenStreetMap (OSM) features based on specific tags.

    This loader uses the `osmnx` library to fetch geographical features from
    OpenStreetMap within a buffered bounding box of the thematic data. It
    automatically handles Coordinate Reference System (CRS) transformations
    between the project CRS and the WGS84 system required by OSM.

    Parameters
    ----------
    osm_tags : dict
        A dictionary of OSM tags used to filter features (e.g., {'highway': True}
        of {'building': 'industrial'}).
    aligner : Aligner
        The aligner object providing the spatial context, target CRS, and logger.
    included_attributes : Iterable[str] | None, optional
        Optional whitelist of OSM attribute names to keep in feature properties.
        - If None (default), all attributes returned by OSMNX are retained.
        - If provided, only these attributes are kept (plus optional directional
          attributes when `include_directional_attributes=True`).
    include_directional_attributes : bool, default False
        If True, ensure common direction-related OSM attributes are included in
        feature properties (e.g. `oneway`, `junction`, lane/turn-direction fields),
        even when an `included_attributes` whitelist is used.

    Attributes
    ----------
    aligner : Aligner
        Reference to the parent aligner object.
    osm_tags : dict
        The tags used for filtering OSM data.
    included_attributes : set[str] | None
        Effective attribute whitelist (or None when all attributes are kept).
    include_directional_attributes : bool
        Flag indicating whether direction-related attributes are always included.
    data_dict_source : dict
        Metadata dictionary tracking the data source ("OSM") and version date.
    versiondate_info : dict
        Dictionary specifying the version date field name and format.
    """

    def __init__(
        self,
        osm_tags: dict,
        aligner: Any,
        *,
        included_attributes: Iterable[str] | None = None,
        include_directional_attributes: bool = False,
    ):
        super().__init__(data_dict={}, data_dict_properties={})
        self.aligner = aligner
        self.osm_tags = osm_tags
        self.included_attributes = (
            set(included_attributes) if included_attributes is not None else None
        )
        self.include_directional_attributes = include_directional_attributes
        self.data_dict_source["source"] = "OSM"
        # This is the overpass API used by default by osmnx. A more
        # abstract reference might be preferable.
        self.data_dict_source["source_url"] = "https://overpass-api.de/api"
        self.versiondate_info = {"name": VERSION_DATE, "format": DATE_FORMAT}

    def load_data(self) -> Any:
        """
        Download and process OSM features.

        The process involves:
        1. Buffering the thematic union to define the search area.
        2. Transforming the search area to WGS84 (EPSG:4326).
        3. Fetching features via `osmnx` based on the provided tags.
        4. Re-projecting the downloaded features back to the Aligner's CRS.
        5. Optionally filtering/stabilizing attribute payload in `properties`.

        Returns
        -------
        Any
            The result of the parent DictLoader's load_data method, containing
            the downloaded geometries and metadata.

        Raises
        ------
        ValueError
            If thematic data has not been loaded into the aligner prior
            to calling this method.

        Notes
        -----
        The search area is expanded using `OSM_MAX_REFERENCE_BUFFER` to ensure
        that reference features partially outside the thematic area are
        fully captured for alignment.

        Attribute handling:
        - Default behavior keeps all attributes returned by OSMNX.
        - Use `included_attributes` to reduce payload and improve downstream
          performance/memory usage.
        - Use `include_directional_attributes=True` for directed-network use
          cases where one-way and lane direction metadata is relevant.


        """
        if not self.aligner.thematic_data:
            raise ValueError("Thematic data not loaded")

        # 1. Define search area with buffer
        geom_union = buffer_pos(
            self.aligner.thematic_data.union, OSM_MAX_REFERENCE_BUFFER
        )
        if geom_union is None or geom_union.is_empty:
            raise ValueError(
                "Reference could not be loaded. Please load thematic data first"
            )

        # 2. Transform to WGS84 for OSMNX query
        gdf = gpd.GeoDataFrame(geometry=[geom_union], crs=self.aligner.crs)
        gdf.to_crs(crs="EPSG:4326", inplace=True)
        geom_union_wgs84 = gdf.geometry.iloc[0]

        # 3. Fetch data from OSM
        osm_data = ox.features_from_bbox(geom_union_wgs84.bounds, tags=self.osm_tags)
        osm_data = osm_data.reset_index(drop=False)
        osm_data = osm_data.rename(columns={"id": "osm_id"})

        # 4. Re-project back to project CRS
        osm_data.to_crs(crs=self.aligner.crs, inplace=True)

        # 5. Populate the dictionary for the DictLoader
        self.data_dict = {}
        self.data_dict_properties = {}
        directional_keys = {
            "oneway",
            "junction",
            "oneway:bicycle",
            "oneway:motor_vehicle",
            "oneway:motorcar",
            "oneway:bus",
            "oneway:hgv",
        }
        if self.include_directional_attributes:
            directional_keys.update(
                {
                    "highway",
                    "lanes",
                    "lanes:forward",
                    "lanes:backward",
                    "turn:lanes",
                    "turn:lanes:forward",
                    "turn:lanes:backward",
                }
            )
        for _, row in osm_data.iterrows():
            osm_id = row["osm_id"]
            self.data_dict[osm_id] = row["geometry"]
            props = row.drop(labels=["geometry"]).to_dict()
            if self.included_attributes is not None:
                props = {k: v for k, v in props.items() if k in self.included_attributes}
            if self.include_directional_attributes:
                for key in directional_keys:
                    if key in row:
                        props[key] = row[key]
            self.data_dict_properties[osm_id] = props

        self.data_dict_source[VERSION_DATE] = datetime.now().strftime(DATE_FORMAT)
        self.aligner.logger.feedback_info(f"OSM downloaded: {self.osm_tags}")

        return super().load_data()
