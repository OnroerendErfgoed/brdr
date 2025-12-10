from typing import Any

import numpy as np
from shapely import GeometryCollection
from shapely import STRtree
from shapely.geometry.base import BaseGeometry

from brdr.geometry_utils import extract_points_lines_from_geometry
from brdr.geometry_utils import safe_unary_union
from brdr.typings import ThematicId


class AlignerFeature:
    def __init__(
        self, brdr_id: str, geometry: BaseGeometry, properties: dict[str, Any]
    ):
        self.brdr_id = brdr_id
        self.geometry = geometry
        self.properties = properties


class AlignerFeatureCollection:
    def __init__(
        self,
        features: dict[ThematicId, AlignerFeature],
        source: dict[str, str] = None,
        is_reference: bool = False,
    ):
        self.source = source or {}
        self.features = features or {}
        self.is_reference = is_reference

        self._union = None
        self._tree = None
        self._elements = None
        self._items = None

    def __getitem__(self, key: ThematicId):
        return self.features[key]


    @property
    def union(self):
        if not self.is_reference:
            raise ValueError("FeatureCollection is not a reference dataset")
        if not self.features:
            raise ValueError("FeatureCollection has no features")
        if self._union is None:
            geoms = [f.geometry for f in self.features.values()]
            self._union = safe_unary_union(geoms)
        return self._union

    @property
    def elements(self):
        if not self.is_reference:
            raise ValueError("FeatureCollection is not a reference dataset")
        if not self.features:
            raise ValueError("FeatureCollection has no features")
        if self._elements is None:
            geoms = [f.geometry for f in self.features.values()]
            self._elements = extract_points_lines_from_geometry(
                GeometryCollection(geoms)
            )
        return self._elements

    @property
    def items(self):
        if not self.is_reference:
            raise ValueError("FeatureCollection is not a reference dataset")
        if not self.features:
            raise ValueError("FeatureCollection has no features")
        if self._items is None:
            self._items = np.array(list(self.features.keys()), dtype=object)
        return self._items

    @property
    def tree(self):
        if not self.is_reference:
            raise ValueError("FeatureCollection is not a reference dataset")
        if not self.features:
            raise ValueError("FeatureCollection has no features")
        if self._tree is None:
            geoms = [f.geometry for f in self.features.values()]
            self._tree = STRtree(geoms)
        return self._tree
