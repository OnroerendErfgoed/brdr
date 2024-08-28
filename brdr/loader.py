import json
from abc import ABC

from brdr.constants import MAX_REFERENCE_BUFFER
from brdr.geometry_utils import get_bbox, buffer_pos,partition
from brdr.grb import get_reference_data_dict_grb_actual, get_collection_grb_fiscal_parcels
import requests as requests
from shapely import make_valid, unary_union
from shapely.geometry import shape
from shapely.geometry.base import BaseGeometry
from brdr.enums import GRBType
from brdr.typings import FeatureCollection


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

    def load_reference_data_grb_actual(self, *, grb_type=GRBType.ADP, partition=partition):
        data_dict, id_property = get_reference_data_dict_grb_actual(
            aligner=self.aligner,
            grb_type=grb_type,
            partition=partition,
        )
        self.aligner.name_reference_id = id_property
        self.aligner.logger.feedback_info(f"GRB downloaded: {grb_type}")
        self.data_dict = data_dict

class GRBFiscalParcelLoader(GeoJsonLoader):
    def __init__(self, year: str, aligner,prt=1000):
        self.aligner = aligner
        self.year = year
        if not self.aligner.dict_thematic:
            raise ValueError("Thematic data not loaded")
        geom_union = buffer_pos(self.aligner._get_thematic_union(),MAX_REFERENCE_BUFFER)
        collection = {}
        if prt < 1:
            collection = get_collection_grb_fiscal_parcels(year=self.year, bbox=get_bbox(geom_union))
        else:
            geoms = partition(geom_union,prt)
            for g in geoms:
                coll = get_collection_grb_fiscal_parcels(year=self.year, bbox=get_bbox(g))
                if collection=={}:
                    collection = dict(coll)
                elif "features" in collection:
                    collection["features"].extend(coll["features"])
        _input = dict(collection)
        self.aligner.logger.feedback_info(f"Adpf downloaded for year: {self.year}")
        super().__init__(_input, "CAPAKEY")
