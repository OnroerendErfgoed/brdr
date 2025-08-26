from datetime import datetime

import geopandas as gpd
import osmnx as ox

from brdr.constants import (
    DATE_FORMAT,
    VERSION_DATE,
    OSM_MAX_REFERENCE_BUFFER,
)
from brdr.geometry_utils import buffer_pos
from brdr.loader import DictLoader


class OSMLoader(DictLoader):
    def __init__(self, osm_tags, aligner):
        super().__init__(data_dict={}, data_dict_properties={})
        self.aligner = aligner
        self.osm_tags = osm_tags
        self.data_dict_source["source"] = "OSM"
        self.versiondate_info = {"name": VERSION_DATE, "format": DATE_FORMAT}

    def load_data(self):
        if not self.aligner.dict_thematic:
            raise ValueError("Thematic data not loaded")
        geom_union = buffer_pos(
            self.aligner.get_thematic_union(), OSM_MAX_REFERENCE_BUFFER
        )

        gdf = gpd.GeoDataFrame(geometry=[geom_union], crs=self.aligner.CRS)
        gdf.to_crs(crs="EPSG:4326", inplace=True)
        geom_union_wgs84 = gdf.geometry.iloc[0]
        osm_data = ox.features_from_bbox(geom_union_wgs84.bounds, tags=self.osm_tags)
        osm_data = osm_data.reset_index(drop=False)
        osm_data = osm_data.rename(columns={"id": "osm_id"})
        osm_data.to_crs(crs=self.aligner.CRS, inplace=True)
        # osm_data.to_file("dataframe.geojson", driver="GeoJSON")
        self.data_dict = {
            row["osm_id"]: row["geometry"] for idx, row in osm_data.iterrows()
        }
        # self.data_dict_properties = {row['osm_id']: row["geometry"] for idx, row in osm_data.iterrows()}
        self.data_dict_source[VERSION_DATE] = datetime.now().strftime(DATE_FORMAT)
        self.aligner.logger.feedback_info(f"OSM downloaded: {self.osm_tags}")
        return super().load_data()
