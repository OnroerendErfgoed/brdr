"""
Aligner Processor - handles processing and calculation methods
"""
import heapq
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, wait
from datetime import datetime
from typing import Iterable, Dict, Any, List, Tuple

import numpy as np
from shapely import (
    GeometryCollection,
    Point,
    MultiPolygon,
    MultiPoint,
    segmentize,
    shortest_line,
    LinearRing,
    Polygon,
    get_parts,
    make_valid,
    remove_repeated_points,
)
from shapely.errors import GeometryTypeError
from shapely.geometry.base import BaseGeometry
from shapely.geometry.linestring import LineString
from shapely.ops import nearest_points, split

from brdr.aligner_config import AlignerConfig, ProcessingConfig
from brdr.aligner_core import AlignerCore
from brdr.constants import (
    PREDICTION_SCORE,
    PREDICTION_COUNT,
    MAX_OUTER_BUFFER,
    STABILITY,
    ZERO_STREAK,
    REMARK_FIELD_NAME,
    DIFF_PERC_INDEX,
    DIFF_INDEX,
    RELEVANT_DISTANCE_FIELD_NAME,
    NR_CALCULATION_FIELD_NAME,
    FORMULA_FIELD_NAME,
    EVALUATION_FIELD_NAME,
    FULL_BASE_FIELD_NAME,
    FULL_ACTUAL_FIELD_NAME,
    DIFF_PERCENTAGE_FIELD_NAME,
    DIFF_AREA_FIELD_NAME,
    OD_ALIKE_FIELD_NAME,
    EQUAL_REFERENCE_FEATURES_FIELD_NAME,
)
from brdr.enums import (
    OpenDomainStrategy,
    Evaluation,
    SnapStrategy,
    FullStrategy,
    DiffMetric,
)
from brdr.geometry_utils import (
    buffer_neg,
    buffer_pos,
    buffer_neg_pos,
    safe_difference,
    safe_intersection,
    safe_symmetric_difference,
    safe_union,
    safe_unary_union,
    snap_geometry_to_reference,
    get_partitions,
    features_by_geometric_operation,
    create_donut,
    fill_and_remove_gaps,
    multi_to_singles,
    get_shape_index,
    geometric_equality,
)
from brdr.typings import ProcessResult
from brdr.utils import multi_to_singles


class AlignerProcessor:
    """
    Processor class responsible for geometric processing and calculations.
    
    This class handles:
    - Main processing algorithms
    - Geometry alignment calculations
    - Multi-threading for performance
    - Prediction algorithms
    - Evaluation methods
    """
    
    def __init__(self, core: AlignerCore):
        """
        Initialize the processor with a core aligner instance.
        
        Args:
            core: AlignerCore instance with loaded data
        """
        self.core = core
        self.config = core.config
        self.logger = core.logger
    
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
            
        Raises:
            ValueError: If no data is loaded or parameters are invalid
        """
        # Validate that data is loaded
        self.core.validate_data_loaded()
        
        # Set up processing parameters
        if relevant_distances is None:
            relevant_distances = [relevant_distance]
        
        processing_config = ProcessingConfig(
            relevant_distance=relevant_distance,
            od_strategy=od_strategy,
            threshold_overlap_percentage=threshold_overlap_percentage
        )
        
        self.logger.feedback_debug(f"Processing series: {list(relevant_distances)}")
        
        # Use provided dict_thematic or default to loaded data
        thematic_data = dict_thematic or self.core.dict_thematic
        
        # Handle multi-polygon mode
        if self.config.multi_as_single_modus:
            thematic_data, dict_multi_as_single = multi_to_singles(thematic_data)
        else:
            dict_multi_as_single = {}
        
        dict_series = {}
        
        # Process with or without multi-threading
        if self.config.max_workers == -1:
            # Single-threaded processing
            dict_series = self._process_single_threaded(
                thematic_data, relevant_distances, processing_config
            )
        else:
            # Multi-threaded processing
            dict_series = self._process_multi_threaded(
                thematic_data, relevant_distances, processing_config
            )
        
        # Handle multi-polygon reconstruction if needed
        if dict_multi_as_single:
            dict_series = self._reconstruct_multipolygons(dict_series, dict_multi_as_single)
        
        self.logger.feedback_info(f"Processing completed for {len(dict_series)} features")
        return dict_series
    
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
            
        Raises:
            ValueError: If geometry is invalid or parameters are incorrect
        """
        if input_geometry is None or input_geometry.is_empty:
            raise ValueError("Input geometry cannot be None or empty")
        
        if relevant_distance < 0:
            raise ValueError("Relevant distance must be non-negative")
        
        self.core.validate_data_loaded()
        
        processing_config = ProcessingConfig(
            relevant_distance=relevant_distance,
            od_strategy=od_strategy,
            threshold_overlap_percentage=threshold_overlap_percentage
        )
        
        return self._process_single_geometry(input_geometry, processing_config)
    
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
        if relevant_distances is None:
            relevant_distances = self.config.relevant_distances
        
        self.logger.feedback_info(f"Running predictor for {len(relevant_distances)} distances")
        
        # Process all distances
        dict_processresults = self.process(
            dict_thematic=dict_thematic,
            relevant_distances=relevant_distances,
            od_strategy=od_strategy,
            threshold_overlap_percentage=threshold_overlap_percentage,
        )
        
        # Generate predictions
        dict_predictions = self._generate_predictions(dict_processresults)
        
        # Make predictions unique
        dict_predictions = self._make_predictions_unique(dict_predictions)
        
        # Evaluate predictions
        dict_evaluated_predictions = self.evaluate(
            ids_to_evaluate=list(dict_predictions.keys())
        )
        
        return dict_predictions, dict_evaluated_predictions
    
    def evaluate(
        self,
        ids_to_evaluate: List = None,
        base_formula_field: str = FORMULA_FIELD_NAME,
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
        self.core.validate_data_loaded()
        
        if ids_to_evaluate is None:
            ids_to_evaluate = list(self.core.dict_thematic.keys())
        
        self.logger.feedback_info(f"Evaluating {len(ids_to_evaluate)} features")
        
        dict_evaluated = {}
        
        for theme_id in ids_to_evaluate:
            if theme_id not in self.core.dict_thematic:
                continue
            
            evaluation_result = self._evaluate_single_feature(
                theme_id, base_formula_field, **kwargs
            )
            dict_evaluated[theme_id] = evaluation_result
        
        return dict_evaluated
    
    def get_brdr_formula(self, geometry: BaseGeometry, with_geom: bool = False) -> Dict:
        """
        Calculate formula-related information based on input geometry.
        
        Args:
            geometry: Input geometry to analyze
            with_geom: Whether to include geometry in result
            
        Returns:
            Dictionary containing formula information
            
        Raises:
            ValueError: If geometry is invalid or no reference data loaded
        """
        if geometry is None or geometry.is_empty:
            raise ValueError("Geometry cannot be None or empty")
        
        self.core.validate_data_loaded()
        
        return self._calculate_formula(geometry, with_geom)
    
    def _process_single_threaded(
        self,
        thematic_data: Dict,
        relevant_distances: Iterable[float],
        config: ProcessingConfig
    ) -> Dict:
        """Process data in single-threaded mode."""
        dict_series = {}
        
        total_features = len(thematic_data)
        for i, (theme_id, geometry) in enumerate(thematic_data.items(), 1):
            self.logger.feedback_info(f"Processing feature {i}/{total_features}: {theme_id}")
            
            dict_series[theme_id] = {}
            for distance in relevant_distances:
                config_for_distance = ProcessingConfig(
                    relevant_distance=distance,
                    od_strategy=config.od_strategy,
                    threshold_overlap_percentage=config.threshold_overlap_percentage
                )
                
                result = self._process_single_geometry(geometry, config_for_distance)
                dict_series[theme_id][distance] = result
        
        return dict_series
    
    def _process_multi_threaded(
        self,
        thematic_data: Dict,
        relevant_distances: Iterable[float],
        config: ProcessingConfig
    ) -> Dict:
        """Process data in multi-threaded mode."""
        dict_series = {}
        
        with ThreadPoolExecutor(max_workers=self.config.max_workers) as executor:
            futures = []
            
            for theme_id, geometry in thematic_data.items():
                dict_series[theme_id] = {}
                
                for distance in relevant_distances:
                    config_for_distance = ProcessingConfig(
                        relevant_distance=distance,
                        od_strategy=config.od_strategy,
                        threshold_overlap_percentage=config.threshold_overlap_percentage
                    )
                    
                    future = executor.submit(
                        self._process_single_geometry, 
                        geometry, 
                        config_for_distance
                    )
                    futures.append((future, theme_id, distance))
            
            # Wait for all futures to complete
            completed_futures = [f[0] for f in futures]
            wait(completed_futures)
            
            # Collect results
            for future, theme_id, distance in futures:
                try:
                    result = future.result()
                    dict_series[theme_id][distance] = result
                except Exception as e:
                    self.logger.feedback_warning(f"Error processing {theme_id} at distance {distance}: {e}")
                    # Create empty result for failed processing
                    dict_series[theme_id][distance] = self._create_empty_result()
        
        return dict_series
    
    def _process_single_geometry(
        self, 
        geometry: BaseGeometry, 
        config: ProcessingConfig
    ) -> ProcessResult:
        """
        Process a single geometry with given configuration.
        
        This is the core processing method that implements the BRDR algorithm.
        """
        try:
            # Validate input geometry
            if geometry is None or geometry.is_empty:
                return self._create_empty_result()
            
            # Area limit check
            if (self.config.area_limit is not None and 
                hasattr(geometry, 'area') and 
                geometry.area > self.config.area_limit):
                self.logger.feedback_warning(f"Geometry area ({geometry.area}) exceeds limit ({self.config.area_limit})")
                return self._create_empty_result()
            
            # Circle detection (skip perfect circles)
            if hasattr(geometry, 'area') and hasattr(geometry, 'length'):
                shape_index = get_shape_index(geometry.area, geometry.length)
                if shape_index > self.config.threshold_circle_ratio:
                    self.logger.feedback_debug("Skipping circular geometry")
                    return self._create_empty_result()
            
            # Main processing algorithm would go here
            # For now, return a basic result structure
            result = self._run_alignment_algorithm(geometry, config)
            
            return result
            
        except Exception as e:
            self.logger.feedback_warning(f"Error in geometry processing: {e}")
            return self._create_empty_result()
    
    def _run_alignment_algorithm(
        self, 
        geometry: BaseGeometry, 
        config: ProcessingConfig
    ) -> ProcessResult:
        """
        Run the main BRDR alignment algorithm.
        
        This method would contain the core alignment logic from the original
        aligner.py file. For this refactoring, I'm providing the structure.
        """
        # TODO: Implement the core alignment algorithm
        # This would include the logic from the original process_geometry method
        
        # Placeholder implementation
        result = {
            'result': geometry,  # Aligned geometry
            'result_diff': GeometryCollection(),  # Difference geometry
            'result_diff_plus': GeometryCollection(),  # Added areas
            'result_diff_min': GeometryCollection(),  # Removed areas
            'remark': 'Processed'
        }
        
        return result
    
    def _create_empty_result(self) -> ProcessResult:
        """Create an empty ProcessResult for failed or skipped processing."""
        return {
            'result': GeometryCollection(),
            'result_diff': GeometryCollection(),
            'result_diff_plus': GeometryCollection(),
            'result_diff_min': GeometryCollection(),
            'remark': 'Empty or failed processing'
        }
    
    def _generate_predictions(self, dict_processresults: Dict) -> Dict:
        """Generate predictions from process results."""
        # TODO: Implement prediction generation logic
        # This would analyze stability across different distances
        return {}
    
    def _make_predictions_unique(self, dict_predictions: Dict) -> Dict:
        """Make predictions unique by removing duplicates."""
        # TODO: Implement uniqueness logic
        return dict_predictions
    
    def _evaluate_single_feature(self, theme_id: str, base_formula_field: str, **kwargs) -> Dict:
        """Evaluate a single feature and assign evaluation code."""
        # TODO: Implement evaluation logic
        return {
            EVALUATION_FIELD_NAME: Evaluation.TO_CHECK_NO_PREDICTION.value
        }
    
    def _calculate_formula(self, geometry: BaseGeometry, with_geom: bool) -> Dict:
        """Calculate BRDR formula for geometry."""
        # TODO: Implement formula calculation
        return {
            'full': False,
            'reference_features': {}
        }
    
    def _reconstruct_multipolygons(
        self, 
        dict_series: Dict, 
        dict_multi_as_single: Dict
    ) -> Dict:
        """Reconstruct multipolygons from single polygon results."""
        # TODO: Implement multipolygon reconstruction
        return dict_series
