# define a typeddict thematic_data with keys name: str and geom: geometry
from typing import Dict
from typing import List
from typing import TypedDict

from shapely.geometry.base import BaseGeometry


class GeoJSONGeometry(TypedDict):
    type: str
    coordinates: List


class Crs(TypedDict):
    type: str
    properties: dict


class Feature(TypedDict):
    type: str
    geometry: GeoJSONGeometry
    properties: dict


class FeatureCollection(TypedDict, total=False):
    type: str
    name: str
    crs: Crs
    features: List[Feature]
    __extra_items__: Dict[str, str]


class ProcessResult(TypedDict, total=False):
    result: BaseGeometry
    result_diff: BaseGeometry
    result_diff_plus: BaseGeometry
    result_diff_min: BaseGeometry
    result_relevant_intersection: BaseGeometry
    result_relevant_diff: BaseGeometry
