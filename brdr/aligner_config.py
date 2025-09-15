"""
Configuration classes for the BRDR Aligner
"""
from dataclasses import dataclass
from typing import List, Optional
import numpy as np

from brdr.constants import (
    DEFAULT_CRS,
    PARTIAL_SNAP_MAX_SEGMENT_LENGTH,
    PARTIAL_SNAP_STRATEGY,
    PARTIAL_SNAPPING,
    RELEVANT_DISTANCE_DECIMALS,
    SNAP_STRATEGY,
    SNAP_MAX_SEGMENT_LENGTH,
    BUFFER_MULTIPLICATION_FACTOR,
    DIFF_METRIC,
)
from brdr.enums import OpenDomainStrategy, SnapStrategy, DiffMetric


@dataclass
class AlignerConfig:
    """
    Configuration class for Aligner to reduce parameter complexity in constructor.
    
    This class encapsulates all configuration parameters for the Aligner,
    making it easier to manage and validate configuration options.
    """
    
    # Core parameters
    relevant_distance: float = 1.0
    relevant_distances: List[float] = None
    crs: str = DEFAULT_CRS
    
    # Strategy parameters
    od_strategy: OpenDomainStrategy = OpenDomainStrategy.SNAP_ALL_SIDE
    snap_strategy: SnapStrategy = SNAP_STRATEGY
    partial_snap_strategy: SnapStrategy = PARTIAL_SNAP_STRATEGY
    diff_metric: DiffMetric = DIFF_METRIC
    
    # Threshold parameters
    threshold_overlap_percentage: int = 50
    threshold_exclusion_area: float = 0
    threshold_exclusion_percentage: float = 0
    threshold_inclusion_percentage: float = 100
    threshold_circle_ratio: float = 0.98
    
    # Processing parameters
    multi_as_single_modus: bool = True
    preserve_topology: bool = False
    partial_snapping: bool = PARTIAL_SNAPPING
    
    # Distance and buffer parameters
    snap_max_segment_length: float = SNAP_MAX_SEGMENT_LENGTH
    partial_snap_max_segment_length: float = PARTIAL_SNAP_MAX_SEGMENT_LENGTH
    buffer_multiplication_factor: float = BUFFER_MULTIPLICATION_FACTOR
    correction_distance: float = 0.01
    mitre_limit: int = 10
    
    # Performance parameters
    area_limit: Optional[float] = None
    max_workers: Optional[int] = None
    
    # UI/Feedback parameters
    feedback: Optional[object] = None
    
    def __post_init__(self):
        """Post-initialization validation and setup."""
        if self.relevant_distances is None:
            self.relevant_distances = [
                round(k, RELEVANT_DISTANCE_DECIMALS)
                for k in np.arange(0, 310, 10, dtype=int) / 100
            ]
        
        self._validate_config()
    
    def _validate_config(self):
        """Validate configuration parameters."""
        if self.relevant_distance < 0:
            raise ValueError("relevant_distance must be non-negative")
        
        if not all(d >= 0 for d in self.relevant_distances):
            raise ValueError("All relevant_distances must be non-negative")
        
        if not (0 <= self.threshold_overlap_percentage <= 100):
            raise ValueError("threshold_overlap_percentage must be between 0 and 100")
        
        if self.threshold_exclusion_area < 0:
            raise ValueError("threshold_exclusion_area must be non-negative")
        
        if not (0 <= self.threshold_exclusion_percentage <= 100):
            raise ValueError("threshold_exclusion_percentage must be between 0 and 100")
        
        if not (0 <= self.threshold_inclusion_percentage <= 100):
            raise ValueError("threshold_inclusion_percentage must be between 0 and 100")
        
        if not (0 <= self.threshold_circle_ratio <= 1):
            raise ValueError("threshold_circle_ratio must be between 0 and 1")
        
        if self.correction_distance < 0:
            raise ValueError("correction_distance must be non-negative")
        
        if self.mitre_limit <= 0:
            raise ValueError("mitre_limit must be positive")


@dataclass 
class ProcessingConfig:
    """Configuration specific to processing operations."""
    
    relevant_distance: float
    od_strategy: OpenDomainStrategy
    threshold_overlap_percentage: int
    
    @classmethod
    def from_aligner_config(cls, config: AlignerConfig, **overrides) -> 'ProcessingConfig':
        """Create ProcessingConfig from AlignerConfig with optional overrides."""
        return cls(
            relevant_distance=overrides.get('relevant_distance', config.relevant_distance),
            od_strategy=overrides.get('od_strategy', config.od_strategy),
            threshold_overlap_percentage=overrides.get('threshold_overlap_percentage', 
                                                       config.threshold_overlap_percentage)
        )
