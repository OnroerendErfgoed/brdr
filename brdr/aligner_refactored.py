"""
Refactored Aligner class - main interface combining core, processor, and exporter
"""
from typing import Dict, Any, Iterable, List, Tuple, Optional

from brdr.aligner_config import AlignerConfig
from brdr.aligner_core import AlignerCore
from brdr.aligner_processor import AlignerProcessor
from brdr.aligner_exporter import AlignerExporter
from brdr.enums import AlignerResultType, AlignerInputType, OpenDomainStrategy
from brdr.loader import Loader
from brdr.typings import ProcessResult
from shapely.geometry.base import BaseGeometry


class Aligner:
    """
    Main Aligner class that provides the complete BRDR alignment functionality.
    
    This refactored version splits the original monolithic Aligner class into
    three focused components:
    - AlignerCore: Data loading and basic operations
    - AlignerProcessor: Processing and calculation methods  
    - AlignerExporter: Export and serialization functionality
    
    The main Aligner class provides a unified interface while maintaining
    backward compatibility with the original API.
    """
    
    def __init__(
        self,
        *,
        feedback=None,
        relevant_distance=1,
        relevant_distances=None,
        threshold_overlap_percentage=50,
        od_strategy=OpenDomainStrategy.SNAP_ALL_SIDE,
        crs="EPSG:31370",
        multi_as_single_modus=True,
        preserve_topology=False,
        snap_strategy=None,
        snap_max_segment_length=2,
        partial_snapping=False,
        partial_snap_strategy=None,
        partial_snap_max_segment_length=2,
        threshold_exclusion_area=0,
        threshold_exclusion_percentage=0,
        threshold_inclusion_percentage=100,
        buffer_multiplication_factor=1.01,
        threshold_circle_ratio=0.98,
        correction_distance=0.01,
        diff_metric=None,
        mitre_limit=10,
        area_limit=None,
        max_workers=None,
    ):
        """
        Initialize the Aligner with configuration parameters.
        
        This constructor maintains backward compatibility with the original
        Aligner class while using the new configuration system internally.
        
        Args:
            feedback: Feedback object for UI integration (e.g., QGIS)
            relevant_distance: Default relevant distance for processing
            relevant_distances: List of relevant distances for batch processing
            threshold_overlap_percentage: Threshold for reference polygon inclusion
            od_strategy: Strategy for handling open domain areas
            crs: Coordinate reference system
            multi_as_single_modus: Whether to split multipolygons for processing
            preserve_topology: Whether to preserve topological relationships
            snap_strategy: Strategy for snapping operations
            snap_max_segment_length: Maximum segment length for snapping
            partial_snapping: Whether to enable partial snapping
            partial_snap_strategy: Strategy for partial snapping
            partial_snap_max_segment_length: Max segment length for partial snapping
            threshold_exclusion_area: Area threshold for excluding reference polygons
            threshold_exclusion_percentage: Percentage threshold for exclusion
            threshold_inclusion_percentage: Percentage threshold for inclusion
            buffer_multiplication_factor: Factor for buffer calculations
            threshold_circle_ratio: Threshold for detecting circular geometries
            correction_distance: Distance for technical corrections
            diff_metric: Metric for measuring differences
            mitre_limit: Limit for mitre joins in buffering
            area_limit: Maximum area for processing geometries
            max_workers: Number of worker threads for parallel processing
        """
        # Create configuration object
        self.config = AlignerConfig(
            feedback=feedback,
            relevant_distance=relevant_distance,
            relevant_distances=relevant_distances,
            threshold_overlap_percentage=threshold_overlap_percentage,
            od_strategy=od_strategy,
            crs=crs,
            multi_as_single_modus=multi_as_single_modus,
            preserve_topology=preserve_topology,
            snap_strategy=snap_strategy or "prefer_vertices",
            snap_max_segment_length=snap_max_segment_length,
            partial_snapping=partial_snapping,
            partial_snap_strategy=partial_snap_strategy or "prefer_vertices",
            partial_snap_max_segment_length=partial_snap_max_segment_length,
            threshold_exclusion_area=threshold_exclusion_area,
            threshold_exclusion_percentage=threshold_exclusion_percentage,
            threshold_inclusion_percentage=threshold_inclusion_percentage,
            buffer_multiplication_factor=buffer_multiplication_factor,
            threshold_circle_ratio=threshold_circle_ratio,
            correction_distance=correction_distance,
            diff_metric=diff_metric or "changes_area",
            mitre_limit=mitre_limit,
            area_limit=area_limit,
            max_workers=max_workers,
        )
        
        # Initialize components
        self.core = AlignerCore(self.config)
        self.processor = AlignerProcessor(self.core)
        self.exporter = AlignerExporter(self.core)
        
        # Expose logger for backward compatibility
        self.logger = self.core.logger
    
    # Data Loading Methods (delegate to core)
    
    def load_thematic_data(self, loader: Loader) -> None:
        """
        Load thematic features into the aligner.
        
        Args:
            loader: Loader instance containing thematic data
        """
        self.core.load_thematic_data(loader)
    
    def load_reference_data(self, loader: Loader) -> None:
        """
        Load reference features into the aligner and prepare for processing.
        
        Args:
            loader: Loader instance containing reference data
        """
        self.core.load_reference_data(loader)
    
    # Processing Methods (delegate to processor)
    
    def process(
        self,
        dict_thematic: Dict = None,
        relevant_distances: Iterable[float] = None,
        relevant_distance: float = 1,
        od_strategy: OpenDomainStrategy = OpenDomainStrategy.SNAP_ALL_SIDE,
        threshold_overlap_percentage: int = 50,
    ) -> Dict[Any, Dict[float, ProcessResult]]:
        """
        Calculate the resulting dictionaries for thematic data based on relevant distances.
        
        Args:
            dict_thematic: Optional override for thematic data
            relevant_distances: Series of relevant distances to process
            relevant_distance: Single relevant distance if not using series
            od_strategy: Strategy for handling open domain areas
            threshold_overlap_percentage: Threshold for reference polygon inclusion
            
        Returns:
            Dictionary with results for each thematic ID and distance
        """
        return self.processor.process(
            dict_thematic=dict_thematic,
            relevant_distances=relevant_distances,
            relevant_distance=relevant_distance,
            od_strategy=od_strategy,
            threshold_overlap_percentage=threshold_overlap_percentage,
        )
    
    def process_geometry(
        self,
        input_geometry: BaseGeometry,
        relevant_distance: float = 1,
        od_strategy: OpenDomainStrategy = OpenDomainStrategy.SNAP_ALL_SIDE,
        threshold_overlap_percentage: int = 50,
    ) -> ProcessResult:
        """
        Process a single geometry with the alignment algorithm.
        
        Args:
            input_geometry: Geometry to process
            relevant_distance: Distance for alignment operations
            od_strategy: Strategy for open domain handling
            threshold_overlap_percentage: Threshold for reference inclusion
            
        Returns:
            ProcessResult containing aligned geometry and metadata
        """
        return self.processor.process_geometry(
            input_geometry=input_geometry,
            relevant_distance=relevant_distance,
            od_strategy=od_strategy,
            threshold_overlap_percentage=threshold_overlap_percentage,
        )
    
    def predictor(
        self,
        dict_thematic: Dict = None,
        relevant_distances: List[float] = None,
        od_strategy: OpenDomainStrategy = OpenDomainStrategy.SNAP_ALL_SIDE,
        threshold_overlap_percentage: int = 50,
    ) -> Tuple[Dict, Dict]:
        """
        Predict stable geometries across multiple relevant distances.
        
        Args:
            dict_thematic: Optional thematic data override
            relevant_distances: Distances to analyze for predictions
            od_strategy: Open domain strategy
            threshold_overlap_percentage: Threshold for reference inclusion
            
        Returns:
            Tuple of (predictions_dict, evaluated_predictions_dict)
        """
        return self.processor.predictor(
            dict_thematic=dict_thematic,
            relevant_distances=relevant_distances,
            od_strategy=od_strategy,
            threshold_overlap_percentage=threshold_overlap_percentage,
        )
    
    def evaluate(
        self,
        ids_to_evaluate: List = None,
        base_formula_field: str = "brdr_formula",
        **kwargs
    ) -> Dict:
        """
        Evaluate predictions and assign evaluation codes.
        
        Args:
            ids_to_evaluate: List of IDs to evaluate (None for all)
            base_formula_field: Field name for base formula comparison
            **kwargs: Additional evaluation parameters
            
        Returns:
            Dictionary with evaluation results
        """
        return self.processor.evaluate(
            ids_to_evaluate=ids_to_evaluate,
            base_formula_field=base_formula_field,
            **kwargs
        )
    
    def get_brdr_formula(self, geometry: BaseGeometry, with_geom: bool = False) -> Dict:
        """
        Calculate formula-related information based on input geometry.
        
        Args:
            geometry: Input geometry to analyze
            with_geom: Whether to include geometry in result
            
        Returns:
            Dictionary containing formula information
        """
        return self.processor.get_brdr_formula(geometry, with_geom)
    
    # Export Methods (delegate to exporter)
    
    def get_results_as_geojson(
        self,
        dict_processresults: Dict[Any, Dict[float, ProcessResult]] = None,
        resulttype: AlignerResultType = AlignerResultType.PROCESSRESULTS,
        formula: bool = False,
        attributes: bool = True,
    ) -> Dict[str, Any]:
        """
        Export processing results as GeoJSON feature collections.
        
        Args:
            dict_processresults: Processing results to export
            resulttype: Type of results to export
            formula: Whether to include formula information
            attributes: Whether to include geometry attributes
            
        Returns:
            Dictionary of GeoJSON feature collections by result type
        """
        return self.exporter.get_results_as_geojson(
            dict_processresults=dict_processresults,
            resulttype=resulttype,
            formula=formula,
            attributes=attributes,
        )
    
    def save_results(
        self,
        path: str,
        dict_processresults: Dict[Any, Dict[float, ProcessResult]] = None,
        resulttype: AlignerResultType = AlignerResultType.PROCESSRESULTS,
        formula: bool = True,
        attributes: bool = True,
    ) -> None:
        """
        Save processing results to files.
        
        Args:
            path: Directory path to save files
            dict_processresults: Processing results to save
            resulttype: Type of results to save
            formula: Whether to include formula information
            attributes: Whether to include geometry attributes
        """
        self.exporter.save_results(
            path=path,
            dict_processresults=dict_processresults,
            resulttype=resulttype,
            formula=formula,
            attributes=attributes,
        )
    
    def get_input_as_geojson(
        self, 
        inputtype: AlignerInputType = AlignerInputType.REFERENCE
    ):
        """
        Get input data as GeoJSON FeatureCollection.
        
        Args:
            inputtype: Type of input data to export
            
        Returns:
            GeoJSON FeatureCollection of input data
        """
        return self.exporter.get_input_as_geojson(inputtype)
    
    def get_diff_metrics(
        self,
        dict_processresults: Dict[Any, Dict[float, ProcessResult]] = None,
        dict_thematic: Dict = None,
    ) -> Dict:
        """
        Calculate difference metrics for processing results.
        
        Args:
            dict_processresults: Processing results to analyze
            dict_thematic: Original thematic data for comparison
            
        Returns:
            Dictionary containing difference metrics
        """
        return self.exporter.get_diff_metrics(
            dict_processresults=dict_processresults,
            dict_thematic=dict_thematic,
        )
    
    # Utility Methods (delegate to core)
    
    def get_thematic_union(self) -> BaseGeometry:
        """
        Get the union of all thematic geometries.
        
        Returns:
            BaseGeometry: Union of all thematic geometries
        """
        return self.core.get_thematic_union()
    
    # Properties for backward compatibility
    
    @property
    def dict_thematic(self) -> Dict[Any, BaseGeometry]:
        """Access to thematic data dictionary."""
        return self.core.dict_thematic
    
    @property
    def dict_reference(self) -> Dict[Any, BaseGeometry]:
        """Access to reference data dictionary."""
        return self.core.dict_reference
    
    @property
    def dict_thematic_properties(self) -> Dict[Any, Dict]:
        """Access to thematic properties dictionary."""
        return self.core.dict_thematic_properties
    
    @property
    def dict_reference_properties(self) -> Dict[Any, Dict]:
        """Access to reference properties dictionary."""
        return self.core.dict_reference_properties
    
    @property
    def CRS(self) -> str:
        """Get the coordinate reference system."""
        return self.core.CRS
    
    @property
    def relevant_distances(self) -> List[float]:
        """Get/set relevant distances."""
        return self.config.relevant_distances
    
    @relevant_distances.setter
    def relevant_distances(self, value: List[float]) -> None:
        """Set relevant distances."""
        self.config.relevant_distances = value
    
    @property
    def od_strategy(self) -> OpenDomainStrategy:
        """Get/set open domain strategy."""
        return self.config.od_strategy
    
    @od_strategy.setter
    def od_strategy(self, value: OpenDomainStrategy) -> None:
        """Set open domain strategy."""
        self.config.od_strategy = value
    
    @property
    def threshold_overlap_percentage(self) -> int:
        """Get/set threshold overlap percentage."""
        return self.config.threshold_overlap_percentage
    
    @threshold_overlap_percentage.setter
    def threshold_overlap_percentage(self, value: int) -> None:
        """Set threshold overlap percentage."""
        self.config.threshold_overlap_percentage = value
    
    # Additional properties for commonly accessed configuration
    @property
    def correction_distance(self) -> float:
        """Get correction distance."""
        return self.config.correction_distance
    
    @property
    def mitre_limit(self) -> int:
        """Get mitre limit."""
        return self.config.mitre_limit
    
    @property
    def partial_snapping(self) -> bool:
        """Get/set partial snapping flag."""
        return self.config.partial_snapping
    
    @partial_snapping.setter
    def partial_snapping(self, value: bool) -> None:
        """Set partial snapping flag."""
        self.config.partial_snapping = value
    
    def __repr__(self) -> str:
        """String representation of the Aligner."""
        return (f"Aligner(thematic_features={len(self.dict_thematic)}, "
                f"reference_features={len(self.dict_reference)}, "
                f"crs={self.CRS})")
