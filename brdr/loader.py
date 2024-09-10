import json
from abc import ABC

import requests as requests
from shapely import make_valid
from shapely.geometry.base import BaseGeometry

from brdr.typings import FeatureCollection
from brdr.utils import geojson_to_dicts


class Loader(ABC):
    def __init__(self):
        self.data_dict: dict[str, BaseGeometry] = {}
        self.data_dict_properties: dict[str, dict] = {}
        self.data_dict_source: dict[str, str] = {}

    def load_data(self):
        self.data_dict = {x: make_valid(self.data_dict[x]) for x in self.data_dict}
        return self.data_dict, self.data_dict_properties, self.data_dict_source


class DictLoader(Loader):
    def __init__(self, data_dict: dict[str:BaseGeometry]):
        # TODO: add dict_properties & dict_source?
        super().__init__()
        self.data_dict = data_dict

    def load_data(self):
        # self._prepare_reference_data()
        return super().load_data()


class GeoJsonLoader(Loader):
    def __init__(
        self,
        *,
        id_property: str = None,
        _input: FeatureCollection = None,
        # data_dict_properties=None, TODO?
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
        self.data_dict, self.data_dict_properties = geojson_to_dicts(
            self.input, self.id_property
        )
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
