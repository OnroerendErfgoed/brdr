# define a typeddict thematic_data with keys name: str and geom: geometry
from typing import Dict
from typing import List
from typing import TypedDict


class GeoJSONGeometry(TypedDict):
    type: str
    coordinates: List[float] | List[List[float]] | List[List[List[float]]] | List[List[List[List[float]]]]


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
