from typing import Any

import numpy as np
from geojson import FeatureCollection
from pyproj import CRS
from shapely import GeometryCollection
from shapely import STRtree
from shapely.geometry.base import BaseGeometry

from brdr.constants import ID_REFERENCE_FIELD_NAME, ID_THEME_FIELD_NAME
from brdr.geometry_utils import extract_points_lines_from_geometry, from_crs
from brdr.geometry_utils import safe_unary_union
from brdr.typings import ThematicId
from brdr.utils import _feature_from_geom


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
        id_fieldname: str = None,
            crs: CRS = None,
        is_reference: bool = False,
    ):
        self.source = source or {}
        self.features = features or {}
        self.is_reference = is_reference
        self._id_fieldname = id_fieldname or None
        self.crs = crs or None
        self._union = None
        self._tree = None
        self._elements = None
        self._items = None

    def __getitem__(self, key: ThematicId):
        return self.features[key]

    @property
    def id_fieldname(self):
        if self._id_fieldname is None:
            if self.is_reference:
                self._id_fieldname = ID_REFERENCE_FIELD_NAME
            else:
                self._id_fieldname = ID_THEME_FIELD_NAME
        return self._id_fieldname

    @property
    def union(self):
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

    def to_geojson(self,geom_attributes=False):
        """
        get a geojson of the input polygons (thematic or reference-polygons)

        Get a GeoJSON (FeatureCollection) from a dictionary of IDs (keys) and geometries (values).

        Args:
            dictionary (dict): Dictionary of geometries.
            crs (str): Coordinate reference system.
            id_field (str): Field name for the ID.
            prop_dict (dict, optional): Dictionary of properties.
            geom_attributes (bool, optional): Whether to include geometry attributes.

        Returns:
            FeatureCollection: The GeoJSON FeatureCollection.
        """
        features = []
        for key, feat in self.features.items():
            properties = feat.properties
            properties[self.id_fieldname] = key
            features.append(_feature_from_geom(feat.geometry, key, properties, geom_attributes))
        crs_geojson = None
        if self.crs is not None:
            crs_geojson = {"type": "name", "properties": {"name": from_crs(self.crs)}}
        geojson = FeatureCollection(features, crs=crs_geojson)
        return geojson
