from abc import ABC, abstractmethod

from shapely.geometry.base import BaseGeometry

from brdr.enums import SnapStrategy
from brdr.typings import ProcessResult


class BaseProcessor(ABC):
    def __init__(self):
        pass

    @abstractmethod
    def process_geometry(self, geometry: BaseGeometry)-> ProcessResult:
        pass

class SnapProcessor(BaseProcessor):
    def process_geometry(
        self,
        input_geometry: BaseGeometry,
        ref_intersections_geoms,
        relevant_distance,
        snap_strategy,
        snap_max_segment_length,
    ) -> ProcessResult:
        pass

class NetworkProcessor(BaseProcessor):
    def process_geometry(
        self,
        input_geometry: BaseGeometry,
        relevant_distance: float = 1,
        snap_strategy=SnapStrategy.NO_PREFERENCE,
    ) -> ProcessResult:
        pass


class DieussaertProcessor(BaseProcessor):
    def process_geometry(
        self,
        input_geometry,
        input_geometry_outer,
        input_geometry_inner,
        ref_intersections_geoms,
        relevant_distance,
    ):
        pass


    
