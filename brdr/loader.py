import json
from abc import ABC

import requests as requests
from shapely import buffer
from shapely import make_valid
from shapely import unary_union
from shapely.geometry import shape
from shapely.geometry.base import BaseGeometry
from shapely.prepared import prep

from brdr.constants import *
from brdr.enums import GRBType
from brdr.geometry_utils import grid_bounds
from brdr.typings import FeatureCollection
from brdr.utils import get_collection


class Loader(ABC):
    def __init__(self):
        self.data_dict: dict[str, BaseGeometry] = {}

    def load_data(self):
        return self.data_dict


class DictLoader(Loader):
    def __init__(self, data_dict: dict[str:BaseGeometry]):
        super().__init__()
        self.data_dict = data_dict

    def load_data(self):
        # self._prepare_reference_data()
        return super().load_data()


class GeoJsonLoader(Loader):
    def __init__(
        self,
        _input: FeatureCollection,
        id_property: str,
    ):
        super().__init__()
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
        for f in self.input["features"]:
            key = str(f["properties"][self.id_property])
            geom = shape(f["geometry"])
            self.data_dict[key] = make_valid(geom)
        return


class GeoJsonFileLoader(GeoJsonLoader):
    def __init__(self, path_to_file, id_property):
        with open(path_to_file, "r") as f:
            _input = json.load(f)
        super().__init__(_input, id_property)


class GeoJsonUrlLoader(GeoJsonLoader):
    def __init__(self, url, id_property):
        _input = requests.get(url).json()
        super().__init__(_input, id_property)


class GRBActualLoader(Loader):
    def __init__(self, grb_type: GRBType, partition: int, aligner):
        super().__init__()
        self.aligner = aligner
        self.grb_type = grb_type
        self.part = partition

    def load_data(self):
        if not self.aligner.dict_thematic:
            raise ValueError("Thematic data not loaded")

        self.load_reference_data_grb_actual(grb_type=self.grb_type, partition=self.part)
        return super().load_data()

    def load_reference_data_grb_actual(self, *, grb_type=GRBType.ADP, partition=0):
        data_dict, id_property = self.get_reference_data_dict_grb_actual(
            grb_type, partition
        )
        self.aligner.name_reference_id = id_property
        self.aligner.logger.feedback_info(f"GRB downloaded: {grb_type}")

        self.data_dict = data_dict

    def get_reference_data_dict_grb_actual(self, grb_type=GRBType.ADP, partition=0):
        """
        Fetches reference data (administrative plots, buildings, or artwork) from the GRB
        API based on thematic data.

        This function retrieves reference data from the Grootschalig Referentie
        Bestand (GRB) depending on the specified `grb_type` (e.g., administrative
        plots (ADP), buildings (GBG), or artwork (KNW)).
        It uses the bounding boxes of the geometries in the loaded thematic data
        (`self.aligner.dict_thematic`) to filter the relevant reference data
        geographically.

        Args:
            grb_type (GRBType, optional): The type of reference data to retrieve.
                Defaults to GRBType.ADP (administrative plots).
            partition (int, optional): If greater than zero, partitions the bounding box
                of the thematic data into a grid before fetching reference data by
                partition. Defaults to 0 (no partitioning).

        Returns:
            tuple: A tuple containing two elements:

                - dict: A dictionary where keys are reference data identifiers
                  (as defined by `name_reference_id`) and values are GeoJSON  geometry
                  objects representing the reference data.
                - str: The name of the reference data identifier property
                  (e.g., "CAPAKEY" for ADP).

        Raises:
            ValueError: If an unsupported `grb_type` is provided.
        """
        if grb_type == GRBType.ADP:
            url_grb = (
                "https://geo.api.vlaanderen.be/GRB/ogc/features/collections/ADP/items?"
            )
            name_reference_id = "CAPAKEY"
        elif grb_type == "gbg":
            url_grb = (
                "https://geo.api.vlaanderen.be/GRB/ogc/features/collections/GBG/items?"
            )
            name_reference_id = "OIDN"
        elif grb_type == GRBType.KNW:
            url_grb = (
                "https://geo.api.vlaanderen.be/GRB/ogc/features/collections/KNW/items?"
            )
            name_reference_id = "OIDN"
        else:
            self.aligner.logger.feedback_info(
                f"type not implemented: {str(grb_type)} -->No reference-data loaded"
            )
            return

        crs = self.aligner.CRS
        limit = DOWNLOAD_LIMIT
        collection = {}
        bounds_array = []

        # Get the bounds of the thematic_data to get the necessary GRB-data
        for key in self.aligner.dict_thematic:
            # buffer them geometry with x m (default 10)
            buffer_value = self.aligner.relevant_distance + MAX_REFERENCE_BUFFER
            geom = buffer(
                self.aligner.dict_thematic[key],
                buffer_value,
                quad_segs=QUAD_SEGMENTS,
                join_style="mitre",
                mitre_limit=MITRE_LIMIT,
            )
            bounds_array.append(geom)
            if partition < 1:
                bbox = str(geom.bounds).strip("()")
                url_grb_bbox = (
                    url_grb
                    + "f=application%2Fgeo%2Bjson&limit="
                    + str(limit)
                    + "&crs="
                    + crs
                    + "&bbox-crs="
                    + crs
                    + "&bbox="
                    + bbox
                )
                self.aligner.logger.feedback_debug(key + "-->" + str(url_grb_bbox))
                coll = self._get_dict_from_url(url_grb_bbox, name_reference_id, limit)
                collection.update(coll)
            if partition > 0:
                geom = unary_union(bounds_array)
                grid = self.partition(geom, partition)
                for g in grid:
                    bbox = str(g.bounds).strip("()")
                    url_grb_bbox = (
                        url_grb
                        + "f=application%2Fgeo%2Bjson&limit="
                        + str(limit)
                        + "&crs="
                        + crs
                        + "&bbox-crs="
                        + crs
                        + "&bbox="
                        + bbox
                    )
                    self.aligner.logger.feedback_debug(key + "-->" + str(url_grb_bbox))
                    coll = self._get_dict_from_url(
                        url_grb_bbox, name_reference_id, limit
                    )
                    collection.update(coll)

        return collection, name_reference_id

    @staticmethod
    def partition(geom, delta):
        """
        Filters a computed grid of partitions (generated by `_grid_bounds`) based on
        intersection with a geometric object (`geom`).

        Args:
            geom (BaseGeometry): The geometric object to check for intersection
                with partitions.
            delta (float): The distance between partitions (same value used in
                `_grid_bounds`).

        Returns:
            list: A filtered list of Polygon objects representing the partitions
                overlapping the original geometric object.
        """
        prepared_geom = prep(geom)
        partitions = grid_bounds(geom, delta)
        filtered_grid = list(filter(prepared_geom.intersects, partitions))
        return filtered_grid

    def _get_dict_from_url(self, input_url, name_reference_id, limit):
        collection = get_collection(input_url, limit)
        dictionary = {}
        if "features" not in collection or len(collection["features"]) == 0:
            return dictionary
        for f in collection["features"]:
            key = str(f["properties"][name_reference_id])
            geom = shape(f["geometry"])
            if key not in collection:
                dictionary[key] = make_valid(geom)
            self.aligner.logger.feedback_debug(key + "-->" + str(geom))
        return dictionary
