from dataclasses import dataclass
from typing import Optional

from brdr.enums import OpenDomainStrategy, DiffMetric, SnapStrategy


@dataclass
class ProcessorConfig:
    """
    Configuration for the geometric preprocessing of polygons in the brdr pipeline.

    Parameters
    ----------
    multi_as_single_modus : bool, default True
        Modus to handle multipolygons. If True, input multipolygons are split into
        single polygons for processing and merged back afterward. If False,
        multipolygons are processed as a single unit.
    max_outer_buffer: int, default 50
        Value that is used to calculate the boundary of a thematic geometry wherefor the calculation has to be done. (inner part is added)
    od_strategy : OpenDomainStrategy, default OpenDomainStrategy.SNAP_ALL_SIDE
        The strategy used to determine how to handle information outside the
        reference polygons (Open Domain).
    threshold_overlap_percentage : int, default 50
        Threshold (%) to determine if a reference polygon should be included
        when no relevant intersections or differences exist. If set to -1,
        the original border is returned in such cases.
    threshold_exclusion_area : int, default 0
        Area in square meters ($m^2$) below which candidate reference polygons
        are excluded from processing.
    threshold_exclusion_percentage : int, default 0
        Overlap percentage below which candidate reference polygons are excluded.
    threshold_inclusion_percentage : int, default 100
        Percentage threshold above which a reference polygon is automatically
        included in the alignment process.
    buffer_multiplication_factor : float, default 1.01
        Factor used to buffer thematic objects when searching for reference
        borders (buffer = factor * relevant_distance).
    threshold_circle_ratio : float, default 0.98
        Threshold to exclude circular geometries based on the Polsby-Popper
        algorithm (where 1.0 is a perfect circle).
    snap_strategy : SnapStrategy, default SnapStrategy.PREFER_VERTICES
        The primary strategy used for snapping geometry vertices to the reference. When alignment is done by 'SnapGeometryProcessor', This strategy will be applied
    snap_max_segment_length : int, default 2
        The maximum segment length allowed during the snapping process. When alignment is done by 'SnapGeometryProcessor', the input geometry (line, lineair ring,...) will be split up by default in parts of max X meter

    partial_snapping : bool, default False
        Whether to allow snapping of individual segments rather than the
        entire geometry.
    partial_snap_strategy : SnapStrategy, default SnapStrategy.PREFER_VERTICES
        The strategy used specifically for partial snapping operations. When snapping of partial geometries (geom_x) is executed, This strategy will be applied.
    partial_snap_max_segment_length : int, default 2
        The maximum segment length allowed during partial snapping. When real snapping of vertices is used, the input geometry will be split up by default in parts of max X meter
    area_limit : int, optional
        Maximum area for processing. If a polygon exceeds this limit, it may
        be skipped or simplified.
    """

    multi_as_single_modus: bool = True
    max_outer_buffer: int = 50
    od_strategy: OpenDomainStrategy = OpenDomainStrategy.SNAP_ALL_SIDE
    threshold_overlap_percentage: int = 50
    threshold_exclusion_area: int = 0
    threshold_exclusion_percentage: int = 0
    threshold_inclusion_percentage: int = 100
    buffer_multiplication_factor: float = 1.01
    threshold_circle_ratio: float = 0.98
    snap_strategy: SnapStrategy = SnapStrategy.PREFER_VERTICES
    snap_max_segment_length: int = 2
    partial_snapping: bool = False
    partial_snap_strategy: SnapStrategy = SnapStrategy.PREFER_VERTICES
    partial_snap_max_segment_length: int = 2
    area_limit: Optional[int] = None


@dataclass
class AlignerConfig:
    """
    Configuration for the alignment and correction of spatial boundaries.

    Parameters
    ----------
    correction_distance : float, default 0.01
        The maximum distance allowed for moving vertices during the
        alignment/correction phase.
    diff_metric : DiffMetric, default DiffMetric.SYMMETRICAL_AREA_CHANGE
        The metric used to calculate the difference or distance between
        the input geometry and the reference.
    mitre_limit : float, default 10
        The miter ratio limit for buffering. Controls the maximum length of
        join corners to prevent spikes in sharp angles. Corners exceeding
        this ratio are beveled.
    max_workers : int, optional
        The maximum number of parallel worker threads/processes to use.
        If None, the system determines the optimal count.
    log_metadata : bool, default False
        If True, detailed metadata about the alignment process will be added to ProcessResult
    add_observations : bool, default False
        If True, observations of the resulting geometry are added to metadata.
    """

    correction_distance: float = 0.01
    diff_metric: DiffMetric = DiffMetric.SYMMETRICAL_AREA_CHANGE
    mitre_limit: float = 10
    max_workers: int = None
    log_metadata: bool = False
    add_observations: bool = False
