from typing import Any, Dict

import numpy as np
from geojson import FeatureCollection
from pyproj import CRS
from shapely import GeometryCollection
from shapely import STRtree
from shapely.geometry.base import BaseGeometry

from brdr.constants import ID_REFERENCE_FIELD_NAME, ID_THEME_FIELD_NAME
from brdr.geometry_utils import extract_points_lines_from_geometry, from_crs
from brdr.geometry_utils import safe_unary_union
from brdr.typings import InputId
from brdr.utils import _feature_from_geom


class AlignerFeature:
    """
    A container for a single thematic or reference feature within the aligner.

    This class wraps a geometry and its associated metadata, providing a
    standardized structure for the alignment process.

    Attributes
    ----------
    data_id: ThematicId
        Source id for the feature.
    brdr_id : str
        The unique identifier for this feature (ThematicId).
    geometry : BaseGeometry
        The shapely geometry object representing the spatial component.
    properties : dict[str, Any]
        A dictionary containing the original attributes of the feature.

    Notes
    -----
    The `AlignerFeature` is typically stored within an `AlignerFeatureCollection`.
    It serves as the atomic unit that the `BaseProcessor` acts upon.



    Examples
    --------
    >>> from shapely.geometry import Point
    >>> feature = AlignerFeature(
    ...     brdr_id="feat_01",
    ...     geometry=Point(0, 0),
    ...     properties={"type": "boundary"}
    ... )
    """

    def __init__(
        self,
        data_id: InputId,
        brdr_id: str,
        geometry: BaseGeometry,
        properties: dict[str, Any],
    ):
        """
        Initializes an AlignerFeature instance.

        Parameters
        ----------
        data_id: ThematicId
            Initial Id of the inputdata (thematic/reference geometries).
        brdr_id : str
            Unique identifier for the feature.
        geometry : BaseGeometry
            The spatial geometry (Point, Line, Polygon, etc.).
        properties : dict[str, Any]
            The attribute data associated with this feature.
        """
        self.data_id = data_id
        self.brdr_id = brdr_id
        self.geometry = geometry
        self.properties = properties


class AlignerFeatureCollection:
    """
    A collection of AlignerFeature objects with spatial indexing and utility properties.

    This class manages a group of features, providing high-level access to
    spatial operations like unary unions and R-tree indexing for fast
    spatial queries.

    Attributes
    ----------
    features : dict[ThematicId, AlignerFeature]
        A dictionary mapping feature IDs to their corresponding AlignerFeature objects.
    source : dict[str, str]
        Metadata regarding the source of the data.
    is_reference : bool
        Flag indicating if this collection serves as the reference dataset
        (enables spatial indexing).
    crs : CRS, optional
        The Coordinate Reference System associated with the collection.

    Notes
    -----
    When `is_reference` is set to True, the collection lazily initializes a
    spatial index ([shapely.strtree.STRtree][]) upon the first access to the
    `tree` property.



    Examples
    --------
    >>> collection = AlignerFeatureCollection(features=my_feature_dict, is_reference=True)
    >>> # Access the spatial index
    >>> index = collection.tree
    >>> # Get a union of all geometries
    >>> total_area = collection.union.area
    """

    def __init__(
        self,
        features: dict[InputId, AlignerFeature],
        source: dict[str, str] = None,
        id_fieldname: str = None,
        crs: CRS = None,
        is_reference: bool = False,
    ):
        """
        Initializes the AlignerFeatureCollection.

        Parameters
        ----------
        features : dict[ThematicId, AlignerFeature]
            Dictionary of features to include in the collection.
        source : dict[str, str], optional
            Source metadata.
        id_fieldname : str, optional
            The name of the identifier field of the geometry. If None, it defaults to
            standard constants based on `is_reference`.
        crs : CRS, optional
            The coordinate system of the features.
        is_reference : bool, optional
            Whether this collection is a reference dataset. Defaults to False.
        """
        self.source = source or {}
        self.features = features or {}
        self.is_reference = is_reference
        self._id_fieldname = id_fieldname or None
        self.crs = crs or None
        self._union = None
        self._tree = None
        self._elements = None
        self._items = None
        self._reference_lookup = None

    def __getitem__(self, key: InputId):
        return self.features[key]

    @property
    def id_fieldname(self) -> str:
        """
        The field name used as the primary identifier.

        Returns
        -------
        str
            The identifier field name, determined by the collection type
            if not explicitly set.
        """
        if self._id_fieldname is None:
            if self.is_reference:
                self._id_fieldname = ID_REFERENCE_FIELD_NAME
            else:
                self._id_fieldname = ID_THEME_FIELD_NAME
        return self._id_fieldname

    @property
    def union(self) -> BaseGeometry:
        """
        The unary union of all feature geometries in the collection.

        Returns
        -------
        BaseGeometry
            A single geometry representing the combined extent of all features.

        """
        if self._union is None:
            geoms = [f.geometry for f in self.features.values()]
            self._union = safe_unary_union(geoms)
        return self._union

    @property
    def elements(self):
        """
        Extracts individual geometric elements (points and lines) from the collection.

        This property flattens complex geometries into their constituent linear
        and point components, which is often required for specific alignment
        engines that operate on boundaries rather than full polygons.

        Returns
        -------
        GeometryCollection
            A collection containing all extracted points and lines from the features.

        Raises
        ------
        ValueError
            If the collection is not marked as a reference dataset (`is_reference=False`)

        Notes
        -----
        This property is lazily evaluated and cached after the first access.
        """
        if not self.is_reference:
            raise ValueError("FeatureCollection is not a reference dataset")
        if self._elements is None:
            geoms = [f.geometry for f in self.features.values()]
            self._elements = extract_points_lines_from_geometry(
                GeometryCollection(geoms)
            )
        return self._elements

    @property
    def items(self):
        """
        A Numpy array containing all feature identifiers (ReferenceIds).

        This property provides a high-performance bridge to Numpy operations,
        allowing for vectorized indexing and fast lookups of feature keys.

        Returns
        -------
        np.ndarray
            A 1D Numpy array of objects containing the keys of the `features` dictionary.

        Raises
        ------
        ValueError
            If the collection is not marked as a reference dataset (`is_reference=False`)

        Notes
        -----
        The use of `dtype=object` ensures that complex string or numeric IDs are
        preserved while benefiting from Numpy's array access speeds.



        Examples
        --------
        >>> # Get all IDs as a numpy array for vectorized operations
        >>> id_array = collection.items
        >>> print(type(id_array))
        <class 'numpy.ndarray'>
        """
        if not self.is_reference:
            raise ValueError("FeatureCollection is not a reference dataset")
        if self._items is None:
            self._items = np.array(list(self.features.keys()), dtype=object)
        return self._items

    @property
    def tree(self) -> STRtree:
        """
        The spatial R-tree index for the collection.

        Returns
        -------
        STRtree
            An STRtree object for fast spatial queries.

        Raises
        ------
        ValueError
            If the collection is not marked as a reference dataset
        """
        if not self.is_reference:
            raise ValueError("FeatureCollection is not a reference dataset")
        if self._tree is None:
            geoms = [f.geometry for f in self.features.values()]
            self._tree = STRtree(geoms)
        return self._tree

    @property
    def reference_lookup(self) -> Dict[Any, Any]:
        """
        A Dictionary Lookup table with the mapping of brdr_id and data_id, that can be used to derive the ID when creating metadata

        Returns
        -------
        Dict
            A Dictionary Lookup table

        Raises
        ------
        ValueError
            If the collection is not marked as a reference dataset
        """
        if not self.is_reference:
            raise ValueError("FeatureCollection is not a reference dataset")
        if self._reference_lookup is None:
            self._reference_lookup = {
                f.data_id: f.brdr_id for f in self.features.values()
            }
        return self._reference_lookup

    def to_geojson(self, geom_attributes=False):
        """
        Converts the collection into a GeoJSON FeatureCollection format.

        Parameters
        ----------
        geom_attributes : bool, optional
            If True, includes calculated geometric attributes (like area or length)
            in the feature properties. Defaults to False.

        Returns
        -------
        FeatureCollection
            A GeoJSON representation of the data.
        """
        features = []
        for key, feat in self.features.items():
            properties = feat.properties
            properties[self.id_fieldname] = key
            features.append(
                _feature_from_geom(feat.geometry, key, properties, geom_attributes)
            )
        crs_geojson = None
        if self.crs is not None:
            crs_geojson = {"type": "name", "properties": {"name": from_crs(self.crs)}}
        geojson = FeatureCollection(features, crs=crs_geojson)
        return geojson
