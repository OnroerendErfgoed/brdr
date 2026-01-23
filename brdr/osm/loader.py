from datetime import datetime
from typing import Any

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

    Attributes
    ----------
    aligner : Aligner
        Reference to the parent aligner object.
    osm_tags : dict
        The tags used for filtering OSM data.
    data_dict_source : dict
        Metadata dictionary tracking the data source ("OSM") and version date.
    versiondate_info : dict
        Dictionary specifying the version date field name and format.
    """

    def __init__(self, osm_tags: dict, aligner: Any):
        super().__init__(data_dict={}, data_dict_properties={})
        self.aligner = aligner
        self.osm_tags = osm_tags
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
        self.data_dict = {
            row["osm_id"]: row["geometry"] for idx, row in osm_data.iterrows()
        }

        self.data_dict_source[VERSION_DATE] = datetime.now().strftime(DATE_FORMAT)
        self.aligner.logger.feedback_info(f"OSM downloaded: {self.osm_tags}")

        return super().load_data()
