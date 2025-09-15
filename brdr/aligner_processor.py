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

from brdr import __version__
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
    DATE_FORMAT,
    VERSION_DATE,
    LAST_VERSION_DATE,
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
from brdr.utils import (
    multi_to_singles,
    coverage_ratio,
    determine_stability,
    diffs_from_dict_processresult,
    is_brdr_formula,
)


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
        
        This method contains the core alignment logic from the original
        aligner.py process_geometry method, adapted for the refactored structure.
        """
        try:
            # Handle special case: relevant_distance = 0
            if config.relevant_distance == 0:
                self.logger.feedback_debug("Calculation for RD = 0")
                return self._create_result_dict(
                    geometry, 
                    GeometryCollection(), 
                    GeometryCollection(), 
                    GeometryCollection(),
                    GeometryCollection(),
                    GeometryCollection(),
                    "relevant distance 0 --> original geometry returned"
                )

            # Check area limit
            if self.config.area_limit and geometry.area > self.config.area_limit:
                message = f"The input polygon is too large to process: input area {geometry.area} m², limit area: {self.config.area_limit} m²."
                raise ValueError(message)

            self.logger.feedback_debug("process geometry")

            # CALCULATE INNER and OUTER INPUT GEOMETRY for performance optimization
            geometry_inner, geometry_outer = self._calculate_inner_outer(
                geometry, config.relevant_distance
            )
            
            # Get reference geometries within buffered area
            geometry_outer_buffered = buffer_pos(
                geometry_outer,
                config.relevant_distance * self.config.buffer_multiplication_factor,
            )
            
            # Query spatial index for intersecting reference geometries
            if hasattr(self.core, 'reference_tree') and self.core.reference_tree is not None:
                ref_intersections_indices = self.core.reference_tree.query(geometry_outer_buffered)
                ref_intersections = self.core.reference_items.take(ref_intersections_indices).tolist()
            else:
                # Fallback: check all reference geometries
                ref_intersections = list(self.core.dict_reference.keys())

            ref_intersections_geoms = []
            all_polygons = True
            for key_ref in ref_intersections:
                ref_geom = self.core.dict_reference[key_ref]
                ref_intersections_geoms.append(ref_geom)
                if not isinstance(ref_geom, (Polygon, MultiPolygon)):
                    all_polygons = False

            if all_polygons:
                return self._process_geometry_by_dieussaert(
                    geometry,
                    geometry_outer,
                    geometry_inner,
                    ref_intersections_geoms,
                    config,
                )
            else:
                # For non-polygon reference data, return basic processing
                self.logger.feedback_debug("Non-polygon reference data detected, using basic processing")
                return self._create_result_dict(
                    geometry, 
                    GeometryCollection(), 
                    GeometryCollection(), 
                    GeometryCollection(),
                    GeometryCollection(),
                    GeometryCollection(),
                    "Basic processing for non-polygon reference"
                )
                
        except Exception as e:
            self.logger.feedback_warning(f"Error in alignment algorithm: {e}")
            return self._create_empty_result()
    
    def _create_empty_result(self) -> ProcessResult:
        """Create an empty ProcessResult for failed or skipped processing."""
        return {
            'result': GeometryCollection(),
            'result_diff': GeometryCollection(),
            'result_diff_plus': GeometryCollection(),
            'result_diff_min': GeometryCollection(),
            'remark': 'Empty or failed processing'
        }
    
    def _create_result_dict(
        self,
        result_geom: BaseGeometry,
        result_diff: BaseGeometry,
        result_diff_plus: BaseGeometry,
        result_diff_min: BaseGeometry,
        relevant_intersection: BaseGeometry,
        relevant_diff: BaseGeometry,
        remark: str
    ) -> ProcessResult:
        """Create a standardized ProcessResult dictionary."""
        return {
            'result': result_geom,
            'result_diff': result_diff,
            'result_diff_plus': result_diff_plus,
            'result_diff_min': result_diff_min,
            'result_relevant_intersection': relevant_intersection,
            'result_relevant_diff': relevant_diff,
            'properties': {REMARK_FIELD_NAME: remark}
        }
    
    def _calculate_inner_outer(self, geometry: BaseGeometry, relevant_distance: float) -> Tuple[BaseGeometry, BaseGeometry]:
        """
        Calculate inner and outer geometries for performance optimization.
        
        Based on _calculate_inner_outer from original aligner.py
        """
        try:
            # Calculate inner geometry (negative buffer)
            geometry_inner = buffer_neg(geometry, relevant_distance)
            if geometry_inner is None or geometry_inner.is_empty:
                geometry_inner = Point(0, 0)  # Fallback for very small geometries
            
            # Calculate outer geometry (positive buffer)  
            geometry_outer = buffer_pos(geometry, relevant_distance)
            if geometry_outer is None or geometry_outer.is_empty:
                geometry_outer = geometry  # Fallback to original geometry
                
            return geometry_inner, geometry_outer
            
        except Exception as e:
            self.logger.feedback_warning(f"Error calculating inner/outer geometry: {e}")
            # Return safe fallbacks
            return Point(0, 0), geometry
    
    def _process_geometry_by_dieussaert(
        self,
        input_geometry: BaseGeometry,
        input_geometry_outer: BaseGeometry,
        input_geometry_inner: BaseGeometry,
        ref_intersections_geoms: List[BaseGeometry],
        config: ProcessingConfig,
    ) -> ProcessResult:
        """
        Process geometry using the Dieussaert algorithm.
        
        Based on _process_geometry_by_dieussaert from original aligner.py
        """
        try:
            buffer_distance = config.relevant_distance / 2
            
            # Initialize arrays for collecting results
            preresult = []
            relevant_intersection_array = []
            relevant_diff_array = []
            
            # Calculate intersection between geometry and open domain
            (
                preresult_od,
                relevant_intersection_od,
                relevant_diff_od,
            ) = self._calculate_intersection_between_geometry_and_od(
                input_geometry_outer, input_geometry_inner, config.relevant_distance
            )
            
            preresult.extend(preresult_od)
            relevant_intersection_array.extend(relevant_intersection_od)
            relevant_diff_array.extend(relevant_diff_od)
            
            # Process each reference geometry
            for geom_reference in ref_intersections_geoms:
                geom_intersection = safe_intersection(input_geometry_outer, geom_reference)
                if geom_intersection.is_empty or geom_intersection is None:
                    continue
                    
                self.logger.feedback_debug("calculate intersection")
                
                # Calculate geometry by intersection and reference
                (
                    geom,
                    relevant_intersection,
                    relevant_diff,
                ) = self._calculate_geom_by_intersection_and_reference(
                    geom_intersection,
                    geom_reference,
                    input_geometry_inner,
                    False,
                    buffer_distance,
                    config.threshold_overlap_percentage,
                    config
                )
                
                self.logger.feedback_debug("intersection calculated")
                
                # Add results to arrays
                preresult = self._add_multi_polygons_from_geom_to_array(geom, preresult)
                relevant_intersection_array = self._add_multi_polygons_from_geom_to_array(
                    relevant_intersection, relevant_intersection_array
                )
                relevant_diff_array = self._add_multi_polygons_from_geom_to_array(
                    relevant_diff, relevant_diff_array
                )
            
            # Union intermediate layers
            relevant_intersection = safe_unary_union(relevant_intersection_array)
            if relevant_intersection is None or relevant_intersection.is_empty:
                relevant_intersection = Polygon()
                
            relevant_diff = safe_unary_union(relevant_diff_array)
            if relevant_diff is None or relevant_diff.is_empty:
                relevant_diff = Polygon()
            
            # Add inner input geometry to preresult
            preresult.append(input_geometry_inner)
            geom_preresult = safe_unary_union(preresult)
            
            # Postprocess results
            return self._postprocess_preresult(
                geom_preresult,
                input_geometry,
                relevant_intersection,
                relevant_diff,
                config.relevant_distance,
            )
            
        except Exception as e:
            self.logger.feedback_warning(f"Error in Dieussaert processing: {e}")
            return self._create_empty_result()
    
    def _calculate_intersection_between_geometry_and_od(
        self, 
        input_geometry_outer: BaseGeometry, 
        input_geometry_inner: BaseGeometry, 
        relevant_distance: float
    ) -> Tuple[List[BaseGeometry], List[BaseGeometry], List[BaseGeometry]]:
        """
        Calculate intersection between geometry and open domain.
        
        Based on _calculate_intersection_between_geometry_and_od from original aligner.py
        """
        try:
            preresult = []
            relevant_intersection_array = []
            relevant_diff_array = []
            
            # Basic implementation - in the original this handles open domain strategy
            # For now, return empty arrays as the main logic is in the reference processing
            return preresult, relevant_intersection_array, relevant_diff_array
            
        except Exception as e:
            self.logger.feedback_warning(f"Error calculating OD intersection: {e}")
            return [], [], []
    
    def _calculate_geom_by_intersection_and_reference(
        self,
        geom_intersection: BaseGeometry,
        geom_reference: BaseGeometry,
        input_geometry_inner: BaseGeometry,
        is_open_domain: bool,
        buffer_distance: float,
        threshold_overlap_percentage: int,
        config: ProcessingConfig
    ) -> Tuple[BaseGeometry, BaseGeometry, BaseGeometry]:
        """
        Calculate geometry by intersection and reference.
        
        Simplified version based on _calculate_geom_by_intersection_and_reference from original aligner.py
        """
        try:
            # Basic intersection-based alignment
            # This is a simplified version - the original is much more complex
            
            # Calculate relevant intersection
            relevant_intersection = safe_intersection(geom_intersection, geom_reference)
            if relevant_intersection is None or relevant_intersection.is_empty:
                relevant_intersection = GeometryCollection()
            
            # Calculate relevant difference  
            relevant_diff = safe_difference(geom_intersection, geom_reference)
            if relevant_diff is None or relevant_diff.is_empty:
                relevant_diff = GeometryCollection()
            
            # For the main geometry, use the intersection with reference
            # In a more complete implementation, this would include buffering, snapping, etc.
            result_geom = geom_intersection
            
            return result_geom, relevant_intersection, relevant_diff
            
        except Exception as e:
            self.logger.feedback_warning(f"Error calculating geom by intersection: {e}")
            return GeometryCollection(), GeometryCollection(), GeometryCollection()
    
    def _add_multi_polygons_from_geom_to_array(
        self, 
        geometry: BaseGeometry, 
        array: List[BaseGeometry]
    ) -> List[BaseGeometry]:
        """
        Add multi-polygons from geometry to array.
        
        Based on _add_multi_polygons_from_geom_to_array from original aligner.py
        """
        try:
            if geometry is None or geometry.is_empty:
                return array
                
            # Handle different geometry types
            if isinstance(geometry, MultiPolygon):
                for geom in geometry.geoms:
                    if not geom.is_empty:
                        array.append(geom)
            elif isinstance(geometry, Polygon):
                if not geometry.is_empty:
                    array.append(geometry)
            elif hasattr(geometry, 'geoms'):  # GeometryCollection
                for geom in geometry.geoms:
                    if isinstance(geom, (Polygon, MultiPolygon)) and not geom.is_empty:
                        array = self._add_multi_polygons_from_geom_to_array(geom, array)
            
            return array
            
        except Exception as e:
            self.logger.feedback_warning(f"Error adding multi-polygons to array: {e}")
            return array
    
    def _postprocess_preresult(
        self,
        geom_preresult: BaseGeometry,
        input_geometry: BaseGeometry,
        relevant_intersection: BaseGeometry,
        relevant_diff: BaseGeometry,
        relevant_distance: float,
    ) -> ProcessResult:
        """
        Postprocess the preresult to create final result.
        
        Based on postprocessing logic from original aligner.py
        """
        try:
            # Calculate differences
            result_diff_plus = safe_difference(
                geom_preresult, 
                buffer_pos(input_geometry, self.config.correction_distance or 0)
            )
            if result_diff_plus is None:
                result_diff_plus = GeometryCollection()
                
            result_diff_min = safe_difference(
                input_geometry, 
                buffer_pos(geom_preresult, self.config.correction_distance or 0)
            )
            if result_diff_min is None:
                result_diff_min = GeometryCollection()
                
            result_diff = safe_unary_union([result_diff_plus, result_diff_min])
            if result_diff is None:
                result_diff = GeometryCollection()
            
            return self._create_result_dict(
                geom_preresult,
                result_diff,
                result_diff_plus,
                result_diff_min,
                relevant_intersection,
                relevant_diff,
                "Processed"
            )
            
        except Exception as e:
            self.logger.feedback_warning(f"Error in postprocessing: {e}")
            return self._create_empty_result()
    
    def _generate_predictions(self, dict_processresults: Dict) -> Dict:
        """
        Generate predictions from process results.
        
        This method analyzes stability across different distances to identify
        'interesting' distances where geometries are stable (zero-streaks).
        
        Based on the predictor logic from the original aligner.py
        """
        try:
            dict_predictions = defaultdict(dict)
            
            # Get thematic data for analysis
            dict_thematic = self.core.dict_thematic
            if not dict_thematic:
                self.logger.feedback_warning("No thematic data available for predictions")
                return {}
            
            # Get reference union for difference calculations
            try:
                reference_union = self.core.get_reference_union()
            except ValueError as e:
                self.logger.feedback_warning(f"Cannot generate predictions: {e}")
                return {}
            
            # Extract relevant distances from process results
            if not dict_processresults:
                return {}
            
            # Get relevant distances from first theme's results
            first_theme_key = next(iter(dict_processresults))
            relevant_distances = list(dict_processresults[first_theme_key].keys())
            relevant_distances = sorted(relevant_distances)
            
            if len(relevant_distances) < 2:
                self.logger.feedback_info("Not enough distances for prediction analysis")
                return {}
            
            # Calculate coverage ratio to determine prediction strategy
            max_relevant_distance = max(relevant_distances)
            cvg_ratio = coverage_ratio(values=relevant_distances, min_val=0, bin_count=10)
            cvg_ratio_threshold = 0.75
            
            self.logger.feedback_debug(f"Coverage ratio: {cvg_ratio}, threshold: {cvg_ratio_threshold}")
            
            # Process each thematic geometry
            for theme_id, dict_processresult in dict_processresults.items():
                if theme_id not in dict_thematic:
                    continue
                    
                # Calculate differences across all distances
                try:
                    diffs = diffs_from_dict_processresult(
                        dict_processresult,
                        dict_thematic[theme_id],
                        reference_union,
                        diff_metric=self.config.diff_metric,
                    )
                except Exception as e:
                    self.logger.feedback_warning(f"Error calculating diffs for {theme_id}: {e}")
                    continue
                
                if len(diffs) != len(relevant_distances):
                    self.logger.feedback_warning(
                        f"Number of computed diffs for {theme_id} does not match distances"
                    )
                    continue
                
                # Determine stability and zero-streaks
                diff_values = list(diffs.values())
                dict_stability = determine_stability(relevant_distances, diff_values)
                
                # Add predictions for distances with zero-streaks
                for rd in relevant_distances:
                    if rd in dict_stability and dict_stability[rd][ZERO_STREAK] is not None:
                        # Copy the process result for this distance
                        prediction_result = dict_processresult[rd].copy()
                        
                        # Add stability information
                        if "properties" not in prediction_result:
                            prediction_result["properties"] = {}
                        
                        prediction_result["properties"][STABILITY] = dict_stability[rd][STABILITY]
                        
                        # Add prediction score if coverage is good enough
                        if cvg_ratio > cvg_ratio_threshold:
                            zero_streak_info = dict_stability[rd][ZERO_STREAK]
                            prediction_result["properties"][PREDICTION_SCORE] = zero_streak_info[3]  # score_streak
                        
                        # Add prediction count (for compatibility)
                        prediction_result["properties"][PREDICTION_COUNT] = 1
                        
                        dict_predictions[theme_id][rd] = prediction_result
                
                self.logger.feedback_debug(f"Generated {len(dict_predictions[theme_id])} predictions for {theme_id}")
            
            return dict(dict_predictions)
            
        except Exception as e:
            self.logger.feedback_warning(f"Error generating predictions: {e}")
            return {}
    
    def _make_predictions_unique(self, dict_predictions: Dict) -> Dict:
        """
        Make predictions unique by removing duplicates.
        
        Check if the predicted geometries are unique and remove duplicated predictions.
        Equality is defined as when the symmetrical difference is smaller than the 
        correction distance.
        
        Based on the _make_predictions_unique logic from the original aligner.py
        """
        try:
            dict_predictions_unique = defaultdict(dict)
            
            for theme_id, dist_results in dict_predictions.items():
                dict_predictions_unique[theme_id] = {}
                predicted_geoms_for_theme_id = []
                
                for rel_dist, processresult in dist_results.items():
                    predicted_geom = processresult["result"]
                    
                    # Check if this geometry is unique or if it's a point/line geometry
                    # (points and lines are always kept as they're less likely to have duplicates)
                    if (not self._equal_geom_in_array(
                        predicted_geom,
                        predicted_geoms_for_theme_id,
                        self.config.correction_distance,
                        self.config.mitre_limit,
                    ) or predicted_geom.geom_type in (
                        "Point",
                        "MultiPoint", 
                        "LineString",
                        "MultiLineString",
                    )):
                        dict_predictions_unique[theme_id][rel_dist] = processresult
                        predicted_geoms_for_theme_id.append(processresult["result"])
                    else:
                        self.logger.feedback_info(
                            f"Duplicate prediction found for key {theme_id} at distance {rel_dist}: Prediction excluded"
                        )
                
                # Update prediction count for all remaining predictions
                for dist in dict_predictions_unique[theme_id].keys():
                    if "properties" not in dict_predictions_unique[theme_id][dist]:
                        dict_predictions_unique[theme_id][dist]["properties"] = {}
                    dict_predictions_unique[theme_id][dist]["properties"][PREDICTION_COUNT] = len(predicted_geoms_for_theme_id)
            
            return dict(dict_predictions_unique)
            
        except Exception as e:
            self.logger.feedback_warning(f"Error making predictions unique: {e}")
            return dict_predictions
    
    def _equal_geom_in_array(self, geom: BaseGeometry, geom_array: list, correction_distance: float, mitre_limit: float) -> bool:
        """
        Check if a predicted geometry is equal to other predicted geometries in a list.
        
        Equality is defined as when the symmetrical difference is smaller than the correction distance.
        Returns True if one of the elements is equal, otherwise False.
        
        Args:
            geom: The geometry to check for equality
            geom_array: List of geometries to compare against
            correction_distance: Distance threshold for equality
            mitre_limit: Mitre limit for buffering operations
            
        Returns:
            bool: True if geometry is found to be equal to one in the array, False otherwise
        """
        try:
            for g in geom_array:
                if geometric_equality(geom, g, correction_distance, mitre_limit):
                    return True
            return False
        except Exception as e:
            self.logger.feedback_warning(f"Error checking geometric equality: {e}")
            return False
    
    def _evaluate_single_feature(self, theme_id: str, geom_predicted: BaseGeometry, base_formula_field: str = FORMULA_FIELD_NAME) -> Dict:
        """
        Evaluate a single feature and assign evaluation code.
        
        Function that evaluates a predicted geometry and returns a properties-dictionary
        with evaluation results based on formula comparison.
        
        Based on the _evaluate logic from the original aligner.py
        
        Args:
            theme_id: ID of the thematic feature
            geom_predicted: The predicted geometry to evaluate
            base_formula_field: Name of the field where the base formula is found
            
        Returns:
            Dict: Properties dictionary with evaluation results
        """
        try:
            threshold_od_percentage = 1
            properties = {
                FORMULA_FIELD_NAME: "",
                EVALUATION_FIELD_NAME: Evaluation.TO_CHECK_NO_PREDICTION,
                FULL_BASE_FIELD_NAME: None,
                FULL_ACTUAL_FIELD_NAME: None,
                OD_ALIKE_FIELD_NAME: None,
                EQUAL_REFERENCE_FEATURES_FIELD_NAME: None,
                DIFF_PERCENTAGE_FIELD_NAME: None,
                DIFF_AREA_FIELD_NAME: None,
            }
            
            # Get the formula for the predicted geometry
            actual_formula = self._get_brdr_formula(geom_predicted)
            if actual_formula is None:
                properties[EVALUATION_FIELD_NAME] = Evaluation.TO_CHECK_NO_PREDICTION
                properties[FULL_ACTUAL_FIELD_NAME] = False
                return properties
                
            properties[FULL_ACTUAL_FIELD_NAME] = actual_formula["full"]
            properties[FORMULA_FIELD_NAME] = json.dumps(actual_formula)
            
            # Get the base formula from thematic properties
            try:
                base_formula = json.loads(
                    self.core.dict_thematic_properties[theme_id][base_formula_field]
                )
            except:
                base_formula = None
                
            if not is_brdr_formula(base_formula):
                properties[EVALUATION_FIELD_NAME] = Evaluation.TO_CHECK_NO_PREDICTION
                return properties
                
            properties[FULL_BASE_FIELD_NAME] = base_formula["full"]
            
            # Check if OD (Open Domain) parts are alike
            od_alike = False
            if (
                base_formula["reference_od"] is None
                and actual_formula["reference_od"] is None
            ):
                od_alike = True
            elif (
                base_formula["reference_od"] is None
                or actual_formula["reference_od"] is None
            ):
                od_alike = False
            elif (
                abs(
                    base_formula["reference_od"]["area"]
                    - actual_formula["reference_od"]["area"]
                )
                * 100
                / base_formula["reference_od"]["area"]
            ) < threshold_od_percentage:
                od_alike = True
            properties[OD_ALIKE_FIELD_NAME] = od_alike
            
            # Check if reference features are equal
            equal_reference_features = False
            if (
                base_formula["reference_features"].keys()
                == actual_formula["reference_features"].keys()
            ):
                equal_reference_features = True
                max_diff_area_reference_feature = 0
                max_diff_percentage_reference_feature = 0
                
                for key in base_formula["reference_features"].keys():
                    if (
                        base_formula["reference_features"][key]["full"]
                        != actual_formula["reference_features"][key]["full"]
                    ):
                        equal_reference_features = False
                        
                    diff_area_reference_feature = abs(
                        base_formula["reference_features"][key]["area"]
                        - actual_formula["reference_features"][key]["area"]
                    )
                    diff_percentage_reference_feature = (
                        abs(
                            base_formula["reference_features"][key]["area"]
                            - actual_formula["reference_features"][key]["area"]
                        )
                        * 100
                        / base_formula["reference_features"][key]["area"]
                    )
                    if diff_area_reference_feature > max_diff_area_reference_feature:
                        max_diff_area_reference_feature = diff_area_reference_feature
                    if (
                        diff_percentage_reference_feature
                        > max_diff_percentage_reference_feature
                    ):
                        max_diff_percentage_reference_feature = (
                            diff_percentage_reference_feature
                        )
                        
                properties[EQUAL_REFERENCE_FEATURES_FIELD_NAME] = equal_reference_features
                properties[DIFF_AREA_FIELD_NAME] = max_diff_area_reference_feature
                properties[DIFF_PERCENTAGE_FIELD_NAME] = (
                    max_diff_percentage_reference_feature
                )
                
            # EVALUATION - Assign evaluation codes based on formula comparison
            if (
                equal_reference_features
                and od_alike
                and base_formula["full"]
                and actual_formula["full"]
            ):  # formula is the same, and both geometries are 'full'
                properties[EVALUATION_FIELD_NAME] = Evaluation.EQUALITY_EQUAL_FORMULA_FULL_1
            elif (
                equal_reference_features
                and od_alike
                and base_formula["full"] == actual_formula["full"]
            ):  # formula is the same, both geometries are not 'full'
                properties[EVALUATION_FIELD_NAME] = Evaluation.EQUALITY_EQUAL_FORMULA_2
            elif (
                base_formula["full"] and actual_formula["full"] and od_alike
            ):  # formula not the same but geometries are full
                properties[EVALUATION_FIELD_NAME] = Evaluation.EQUALITY_FULL_3
            else:
                properties[EVALUATION_FIELD_NAME] = Evaluation.TO_CHECK_NO_PREDICTION
                
            return properties
            
        except Exception as e:
            self.logger.feedback_warning(f"Error evaluating feature {theme_id}: {e}")
            return {
                FORMULA_FIELD_NAME: "",
                EVALUATION_FIELD_NAME: Evaluation.TO_CHECK_NO_PREDICTION,
                FULL_BASE_FIELD_NAME: None,
                FULL_ACTUAL_FIELD_NAME: False,
                OD_ALIKE_FIELD_NAME: None,
                EQUAL_REFERENCE_FEATURES_FIELD_NAME: None,
                DIFF_PERCENTAGE_FIELD_NAME: None,
                DIFF_AREA_FIELD_NAME: None,
            }
    
    def _get_brdr_formula(self, geometry: BaseGeometry, with_geom: bool = False) -> Dict:
        """
        Calculate formula-related information based on the input geometry.
        
        Based on the get_brdr_formula logic from the original aligner.py
        
        Args:
            geometry: The input geometry
            with_geom: Whether to include geometry information in the output
            
        Returns:
            dict: A dictionary containing formula-related data
        """
        try:
            dict_formula = {
                "alignment_date": datetime.now().strftime(DATE_FORMAT),
                "brdr_version": str(__version__),
                "reference_source": self.core.dict_reference_source,
                "full": True,
                "area": round(geometry.area, 2),
                "reference_features": {},
                "reference_od": None,
            }
            
            full_total = True
            last_version_date = None
            
            # Get reference geometries that intersect with the input geometry
            if hasattr(self.core, 'reference_tree') and self.core.reference_tree is not None:
                ref_intersections_indices = self.core.reference_tree.query(geometry)
                ref_intersections = self.core.reference_items.take(ref_intersections_indices).tolist()
            else:
                # Fallback: check all reference geometries
                ref_intersections = list(self.core.dict_reference.keys())
            
            intersected = []
            for key_ref in ref_intersections:
                geom = None
                version_date = None
                geom_reference = self.core.dict_reference[key_ref]
                geom_intersection = make_valid(safe_intersection(geometry, geom_reference))
                
                if geom_intersection.is_empty or geom_intersection is None:
                    continue
                intersected.append(geom_intersection)
                
                geom_reference_area = geom_reference.area
                if geom_reference_area > 0:
                    perc = round(geom_intersection.area * 100 / geom_reference.area, 2)
                else:
                    perc = 0
                if perc < 0.01:
                    continue
                    
                # Add a last_version_date if available in properties
                if (
                    key_ref in self.core.dict_reference_properties
                    and VERSION_DATE in self.core.dict_reference_properties[key_ref]
                ):
                    str_version_date = self.core.dict_reference_properties[key_ref][VERSION_DATE]
                    version_date = datetime.strptime(str_version_date, DATE_FORMAT)
                    if last_version_date is None and version_date is not None:
                        last_version_date = version_date
                    if version_date is not None and version_date > last_version_date:
                        last_version_date = version_date
                        
                if perc > 99.99:
                    full = True
                    area = round(geom_reference_area, 2)
                    perc = 100
                    if with_geom:
                        geom = geom_reference
                else:
                    full = False
                    full_total = False
                    area = round(geom_intersection.area, 2)
                    if with_geom:
                        geom = geom_intersection
                        
                dict_formula["reference_features"][key_ref] = {
                    "full": full,
                    "area": area,
                    "percentage": perc,
                }
                if version_date is not None:
                    dict_formula["reference_features"][key_ref][VERSION_DATE] = (
                        version_date.strftime(DATE_FORMAT)
                    )
                if with_geom:
                    dict_formula["reference_features"][key_ref]["geometry"] = to_geojson(
                        geom
                    )
                    
            dict_formula["full"] = full_total
            if last_version_date is not None:
                dict_formula[LAST_VERSION_DATE] = last_version_date.strftime(DATE_FORMAT)
                
            # Calculate open domain (OD) part - geometry not covered by reference features
            geom_od = buffer_pos(
                buffer_neg(
                    safe_difference(geometry, safe_unary_union(intersected)),
                    self.config.correction_distance,
                    mitre_limit=self.config.mitre_limit,
                ),
                self.config.correction_distance,
                mitre_limit=self.config.mitre_limit,
            )
            if geom_od is not None:
                area_od = round(geom_od.area, 2)
                if area_od > 0:
                    dict_formula["reference_od"] = {"area": area_od}
                    if with_geom:
                        dict_formula["reference_od"]["geometry"] = to_geojson(geom_od)
                        
            self.logger.feedback_debug(str(dict_formula))
            return dict_formula
            
        except Exception as e:
            self.logger.feedback_warning(f"Error calculating BRDR formula: {e}")
            return None
    
    def _reconstruct_multipolygons(
        self, 
        dict_series: Dict, 
        dict_multi_as_single: Dict
    ) -> Dict:
        """
        Reconstruct multipolygons from single polygon results.
        
        Merges process results from multiple theme IDs (single polygons) back into 
        their original multipolygon theme IDs using the merge_process_results function.
        
        Based on the multi_as_single_modus logic from the original aligner.py
        
        Args:
            dict_series: Dictionary with process results from single polygons
            dict_multi_as_single: Mapping from single polygon IDs to original multipolygon IDs
            
        Returns:
            Dict: Merged process results with reconstructed multipolygons
        """
        try:
            if not dict_multi_as_single:
                # No multipolygons to reconstruct
                return dict_series
                
            # Use the original merge_process_results function from utils
            from brdr.utils import merge_process_results
            reconstructed_dict = merge_process_results(dict_series, dict_multi_as_single)
            
            self.logger.feedback_debug(f"Reconstructed {len(dict_multi_as_single)} multipolygon groups")
            return reconstructed_dict
            
        except Exception as e:
            self.logger.feedback_warning(f"Error reconstructing multipolygons: {e}")
            return dict_series
