from dataclasses import dataclass

from brdr.enums import OpenDomainStrategy
from brdr.enums import SnapStrategy


@dataclass
class ProcessorConfig:
    """
    Attributes:
            multi_as_single_modus (boolean, optional): Modus to handle multipolygons (Default=True):
                True: input-multipolygons will be split-up into single polygons and handled by the algorithm. After executing the algorithm, the results are merged together.
                False: Multipolygons are directly processed by the algorithm
        od_strategy (int, optional): The strategy to determine how to handle
            information outside the reference polygons (Open Domain)
            (default: SNAP_FULL_AREA_ALL_SIDE)
        threshold_overlap_percentage (int, optional): Threshold (%) to determine
            from which overlapping-percentage a reference-polygon has to be included
            when there aren't relevant intersections or relevant differences
            (default 50%).
            When setting this parameter to '-1' the original border for will be returned for cases where nor relevant intersections and relevant differences are found
        threshold_exclusion_percentage (int, optional): Percentage for excluding candidate reference-polygons when overlap(%) is smaller than the threshold(Default=0)
        threshold_exclusion_area (int, optional):Area in m² for excluding candidate reference-polygons when overlap(m²) is smaller than the threshold (Default=0)
        buffer_multiplication_factor (float, optional): Multiplication factor, used to buffer the thematic objects searching for reference borders (buffer= buffer_multiplication_factor*relevant_distance)(Default=1.01)
        threshold_circle_ratio (float, optional): Threshold-value to exclude circles getting processed (perfect circle = 1) based on POLSBY-POPPER algorithm(Default=0.98)
        mitre_limit (int, optional):buffer-parameter - The mitre ratio is the ratio of the distance from the corner to the end of the mitred offset corner.
            When two line segments meet at a sharp angle, a miter join will extend far beyond the original geometry. (and in the extreme case will be infinitely far.) To prevent unreasonable geometry, the mitre limit allows controlling the maximum length of the join corner.
            Corners with a ratio which exceed the limit will be beveled(Default=10)
        area_limit (int, optional): Maximum area for processing. (default 100000)
    """

    multi_as_single_modus: bool = True
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
    #TODO: is area limit een processConfig of iets van de aligner?
    area_limit: int = None
