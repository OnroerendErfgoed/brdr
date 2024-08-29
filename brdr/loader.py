import json
from abc import ABC

from brdr.constants import MAX_REFERENCE_BUFFER
from brdr.geometry_utils import buffer_pos
from brdr.grb import get_collection_grb_fiscal_parcels, get_collection_grb_actual
import requests as requests
from shapely import make_valid, unary_union
from shapely.geometry import shape
from shapely.geometry.base import BaseGeometry
from brdr.enums import GRBType
from brdr.typings import FeatureCollection
from brdr.utils import collection_to_dict


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
        self,*,
        id_property: str = None,
            _input: FeatureCollection = None,
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
        self.data_dict = collection_to_dict(self.input, self.id_property)
        return


class GeoJsonFileLoader(GeoJsonLoader):
    def __init__(self, path_to_file, id_property):
        with open(path_to_file, "r") as f:
            _input = json.load(f)
        super().__init__(_input=_input, id_property=id_property)


class GeoJsonUrlLoader(GeoJsonLoader):
    def __init__(self, url, id_property):
        _input = requests.get(url).json()
        super().__init__(_input=_input, id_property=id_property)

class GRBActualLoader(GeoJsonLoader):
    def __init__(self, grb_type: GRBType,aligner, partition: int=1000):
        self.aligner = aligner
        self.grb_type = grb_type
        self.part = partition
        super().__init__()

    def load_data(self):
        if not self.aligner.dict_thematic:
            raise ValueError("Thematic data not loaded")
        geom_union = buffer_pos(self.aligner._get_thematic_union(), MAX_REFERENCE_BUFFER)
        collection, id_property = get_collection_grb_actual(grb_type=self.grb_type, geometry=geom_union, partition=self.part,crs=self.aligner.CRS)
        self.id_property = id_property
        self.input = dict(collection)
        self.aligner.logger.feedback_info(f"GRB downloaded: {self.grb_type}")
        return super().load_data()

class GRBFiscalParcelLoader(GeoJsonLoader):
    def __init__(self, year: str, aligner,partition=1000):
        self.aligner = aligner
        self.year = year
        self.part = partition
        super().__init__(_input=None, id_property="CAPAKEY")

    def load_data(self):
        if not self.aligner.dict_thematic:
            raise ValueError("Thematic data not loaded")
        geom_union = buffer_pos(self.aligner._get_thematic_union(),MAX_REFERENCE_BUFFER)
        collection = get_collection_grb_fiscal_parcels(year=self.year, geometry=geom_union, partition=self.part, crs=self.aligner.CRS)
        self.input = dict(collection)
        self.aligner.logger.feedback_info(f"Adpf downloaded for year: {self.year}")
        return super().load_data()
