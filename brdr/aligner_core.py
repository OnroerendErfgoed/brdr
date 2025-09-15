"""
Core Aligner class - handles data loading and basic alignment operations
"""
import numpy as np
from shapely import STRtree
from shapely.geometry.base import BaseGeometry
from typing import Dict, Any, Optional

from brdr.aligner_config import AlignerConfig
from brdr.geometry_utils import safe_unary_union
from brdr.loader import Loader
from brdr.logger import Logger
from brdr.enums import AlignerInputType
from brdr.utils import geojson_from_dict


class AlignerCore:
    """
    Core Aligner class responsible for data loading and basic operations.
    
    This class handles:
    - Loading thematic and reference data
    - Managing spatial indexes
    - Basic geometric operations
    - Data validation
    """
    
    def __init__(self, config: AlignerConfig):
        """
        Initialize the AlignerCore with configuration.
        
        Args:
            config: AlignerConfig object containing all configuration parameters
        """
        self.config = config
        self.logger = Logger(config.feedback)
        
        # Data storage
        self.dict_thematic: Dict[Any, BaseGeometry] = {}
        self.dict_thematic_properties: Dict[Any, Dict] = {}
        self.dict_thematic_source: Dict[Any, str] = {}
        
        self.dict_reference: Dict[Any, BaseGeometry] = {}
        self.dict_reference_properties: Dict[Any, Dict] = {}
        self.dict_reference_source: Dict[Any, str] = {}
        
        # Spatial indexes and cached unions
        self.reference_tree: Optional[STRtree] = None
        self.reference_items: Optional[np.ndarray] = None
        self.thematic_union: Optional[BaseGeometry] = None
        self.reference_union: Optional[BaseGeometry] = None
        
        # Processing state
        self.relevant_distances = config.relevant_distances
        self.od_strategy = config.od_strategy
        self.threshold_overlap_percentage = config.threshold_overlap_percentage
        
        self.logger.feedback_info(f"AlignerCore initialized with CRS: {config.crs}")
    
    def load_thematic_data(self, loader: Loader) -> None:
        """
        Load thematic features into the aligner.
        
        Args:
            loader: Loader instance containing thematic data
            
        Raises:
            ValueError: If loader is None or contains no data
        """
        if loader is None:
            raise ValueError("Loader cannot be None")
        
        self.logger.feedback_info("Loading thematic data...")
        
        data_dict, data_dict_properties, data_dict_source = loader.load_data()
        
        if not data_dict:
            raise ValueError("Loader contains no thematic data")
        
        self.dict_thematic = data_dict
        self.dict_thematic_properties = data_dict_properties
        self.dict_thematic_source = data_dict_source
        
        # Clear cached union since thematic data changed
        self.thematic_union = None
        
        self.logger.feedback_info(f"Loaded {len(self.dict_thematic)} thematic features")
    
    def load_reference_data(self, loader: Loader) -> None:
        """
        Load reference features into the aligner and prepare for processing.
        
        Args:
            loader: Loader instance containing reference data
            
        Raises:
            ValueError: If loader is None or contains no data
        """
        if loader is None:
            raise ValueError("Loader cannot be None")
        
        self.logger.feedback_info("Loading reference data...")
        
        data_dict, data_dict_properties, data_dict_source = loader.load_data()
        
        if not data_dict:
            raise ValueError("Loader contains no reference data")
        
        self.dict_reference = data_dict
        self.dict_reference_properties = data_dict_properties
        self.dict_reference_source = data_dict_source
        
        # Prepare reference data for processing
        self._prepare_reference_data()
        
        self.logger.feedback_info(f"Loaded {len(self.dict_reference)} reference features")
    
    def get_thematic_union(self) -> BaseGeometry:
        """
        Get the union of all thematic geometries.
        
        Returns:
            BaseGeometry: Union of all thematic geometries
            
        Raises:
            ValueError: If no thematic data is loaded
        """
        if not self.dict_thematic:
            raise ValueError("No thematic data loaded")
        
        if self.thematic_union is None:
            self.logger.feedback_info("Calculating thematic union...")
            self.thematic_union = safe_unary_union(list(self.dict_thematic.values()))
        
        return self.thematic_union
    
    def get_reference_union(self) -> BaseGeometry:
        """
        Get the union of all reference geometries.
        
        Returns:
            BaseGeometry: Union of all reference geometries
            
        Raises:
            ValueError: If no reference data is loaded
        """
        if not self.dict_reference:
            raise ValueError("No reference data loaded")
        
        if self.reference_union is None:
            self.logger.feedback_info("Calculating reference union...")
            self.reference_union = safe_unary_union(list(self.dict_reference.values()))
        
        return self.reference_union
    
    def get_reference_elements(self) -> tuple:
        """
        Get reference elements for spatial queries.
        
        Returns:
            tuple: (reference_tree, reference_items) for spatial operations
            
        Raises:
            ValueError: If no reference data is loaded
        """
        if not self.dict_reference:
            raise ValueError("No reference data loaded")
        
        if self.reference_tree is None or self.reference_items is None:
            self._prepare_reference_data()
        
        return self.reference_tree, self.reference_items
    
    def get_input_as_geojson(self, inputtype: AlignerInputType = AlignerInputType.REFERENCE):
        """
        Get input data as GeoJSON FeatureCollection.
        
        Args:
            inputtype: Type of input data to export (REFERENCE or THEMATIC)
            
        Returns:
            FeatureCollection: GeoJSON representation of the input data
            
        Raises:
            ValueError: If requested input type has no data loaded
        """
        if inputtype == AlignerInputType.REFERENCE:
            if not self.dict_reference:
                raise ValueError("No reference data loaded")
            return geojson_from_dict(
                self.dict_reference, 
                self.config.crs, 
                "ref_id", 
                self.dict_reference_properties
            )
        
        elif inputtype == AlignerInputType.THEMATIC:
            if not self.dict_thematic:
                raise ValueError("No thematic data loaded")
            return geojson_from_dict(
                self.dict_thematic, 
                self.config.crs, 
                "theme_id", 
                self.dict_thematic_properties
            )
        
        else:
            raise ValueError(f"Unknown input type: {inputtype}")
    
    def validate_data_loaded(self) -> None:
        """
        Validate that both thematic and reference data are loaded.
        
        Raises:
            ValueError: If either thematic or reference data is not loaded
        """
        if not self.dict_thematic:
            raise ValueError("No thematic data loaded. Call load_thematic_data() first.")
        
        if not self.dict_reference:
            raise ValueError("No reference data loaded. Call load_reference_data() first.")
    
    def _prepare_reference_data(self) -> None:
        """
        Prepare reference data for spatial queries and analysis.
        
        Creates spatial indexes and clears cached unions.
        """
        if not self.dict_reference:
            return
        
        self.logger.feedback_info(f"Preparing reference data ({len(self.dict_reference)} features)...")
        
        # Create spatial index for performance optimization
        self.reference_tree = STRtree(list(self.dict_reference.values()))
        self.reference_items = np.array(list(self.dict_reference.keys()), dtype=object)
        
        # Clear cached reference union so it will be recalculated when needed
        self.reference_union = None
        
        self.logger.feedback_info("Reference data prepared successfully")
    
    @property
    def CRS(self) -> str:
        """Get the coordinate reference system."""
        return self.config.crs
    
    def __repr__(self) -> str:
        """String representation of AlignerCore."""
        return (f"AlignerCore(thematic_features={len(self.dict_thematic)}, "
                f"reference_features={len(self.dict_reference)}, "
                f"crs={self.config.crs})")

