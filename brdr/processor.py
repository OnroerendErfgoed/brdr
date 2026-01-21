from abc import ABC
from abc import abstractmethod
from typing import List, Any

from shapely import GeometryCollection
from shapely import LinearRing
from shapely import MultiLineString
from shapely import MultiPoint
from shapely import MultiPolygon
from shapely import Point
from shapely import Polygon
from shapely import get_parts
from shapely import make_valid
from shapely import remove_repeated_points
from shapely.errors import GeometryTypeError
from shapely.geometry.base import BaseGeometry
from shapely.geometry.linestring import LineString
from shapely.ops import nearest_points
from shapely.ops import split

from brdr.configs import ProcessorConfig
from brdr.constants import REMARK_FIELD_NAME
from brdr.enums import OpenDomainStrategy
from brdr.enums import ProcessRemark
from brdr.enums import ProcessorID
from brdr.enums import SnapStrategy
from brdr.feature_data import AlignerFeatureCollection
from brdr.geometry_utils import (
    buffer_neg,
    build_custom_network,
)
from brdr.geometry_utils import buffer_neg_pos
from brdr.geometry_utils import buffer_pos
from brdr.geometry_utils import fill_and_remove_gaps
from brdr.geometry_utils import find_best_path_in_network
from brdr.geometry_utils import geometric_equality
from brdr.geometry_utils import get_coords_from_geometry
from brdr.geometry_utils import get_shape_index
from brdr.geometry_utils import safe_difference
from brdr.geometry_utils import safe_intersection
from brdr.geometry_utils import safe_symmetric_difference
from brdr.geometry_utils import safe_unary_union
from brdr.geometry_utils import safe_union
from brdr.geometry_utils import snap_geometry_to_reference
from brdr.geometry_utils import to_multi
from brdr.logger import Logger
from brdr.topo_utils import _dissolve_topo, _generate_topo, _topojson_id_to_arcs
from brdr.typings import ProcessResult, InputId
from brdr.utils import (
    get_relevant_polygons_from_geom,
    build_reverse_index_wkb,
    flatten_iter,
)
from brdr.utils import union_process_result


class BaseProcessor(ABC):
    """
    Abstract base class for geometric alignment processors.

    The processor is responsible for the core alignment logic between thematic
    geometries and reference data. It provides a standardized pipeline for
    processing, cleaning, and validating resulting geometries.

    Attributes
    ----------
    config : ProcessorConfig
        Configuration object containing thresholds and strategy settings.
    feedback : Any, optional
        Feedback mechanism for reporting progress or logs.
    logger : Logger
        Internal logger for debugging and user feedback.
    """

    def __init__(self, config: ProcessorConfig, feedback: Any = None):
        """
        Initializes the BaseProcessor.

        Parameters
        ----------
        config : ProcessorConfig
            Configuration settings for the alignment engine.
        feedback : Any, optional
            Optional feedback handler for the processing status.
        """
        self.feedback = feedback
        self.logger = Logger(feedback)
        self.config = config

    @abstractmethod
    def process(
        self,
        *,
        correction_distance: float,
        reference_data: AlignerFeatureCollection,
        input_geometry: BaseGeometry,
        mitre_limit: float,
        relevant_distance: float,
        thematic_data: AlignerFeatureCollection,
        id_thematic: InputId,
    ) -> ProcessResult:
        """
        Abstract method to align a single geometry to the reference layer.

        Parameters
        ----------
        correction_distance : float
            Small buffer distance used to filter out geometric noise.
        reference_data : AlignerFeatureCollection
            The collection of reference geometries (the 'target').
        input_geometry : BaseGeometry
            The thematic geometry to be aligned.
        mitre_limit : float
            The mitre limit used for buffering operations.
        relevant_distance : float
            The maximum distance within which alignment should occur.
        thematic_data : AlignerFeatureCollection
            The full thematic collection for context.
        id_thematic : Any
            The unique identifier of the feature being processed.

        Returns
        -------
        ProcessResult
            A dictionary containing the resulting geometry and difference metrics.
        """
        pass

    def check_area_limit(self, input_geometry: BaseGeometry):
        if self.config.area_limit and input_geometry.area > self.config.area_limit:
            message = f"The input polygon is too large to process: input area {str(input_geometry.area)} m², limit area: {str(self.config.area_limit)} m²."
            raise ValueError(message)

    def _postprocess_preresult(
        self,
        geom_preresult: BaseGeometry,
        geom_thematic: BaseGeometry,
        relevant_intersection: BaseGeometry,
        relevant_diff: BaseGeometry,
        relevant_distance: float,
        reference_union: BaseGeometry,
        mitre_limit: float,
        correction_distance: float,
    ) -> ProcessResult:
        """
        Refines the initial alignment result through a cleaning pipeline.

        This method applies a series of geometric corrections to ensure the output
        is valid, free of slivers, and respects the specified tolerances.

        Parameters
        ----------
        geom_preresult : BaseGeometry
            The raw geometry result from the alignment engine.
        geom_thematic : BaseGeometry
            The original input geometry.
        relevant_intersection : BaseGeometry
            The intersection part relevant to the alignment.
        relevant_diff : BaseGeometry
            The difference part relevant to the alignment.
        relevant_distance : float
            The maximum search distance used.
        reference_union : BaseGeometry
            The union of all reference features for clipping/context.
        mitre_limit : float
            Limit for sharp corners in buffering.
        correction_distance : float
            Small distance for noise reduction.

        Returns
        -------
        ProcessResult
            A dictionary containing the cleaned 'result' and its differences
            compared to the original input.

        Notes
        -----
        The post-processing pipeline performs the following steps:
        1. **Empty Check**: Returns empty collection if result is null.
        2. **Type Validation**: Ensures geometry type remains consistent.
        3. **Stability Checks**: Reverts to original if it's a circle or unchanged.
        4. **Noise Filtering**: Removes slivers and gaps smaller than `correction_distance`.
        5. **Hole Filling**: Cleans up internal 'donuts' within the geometry.
        6. **Diff Calculation**: Generates `result_diff_plus` and `result_diff_min`.



        ```{mermaid}
        graph TD
            Pre[Pre-result] --> Valid{Is Valid?}
            Valid -- No --> Empty[Return Empty]
            Valid -- Yes --> Circle{Is Circle?}
            Circle -- Yes --> Orig[Return Original]
            Circle -- No --> Clean[Sliver & Gap Removal]
            Clean --> Diff[Calculate Differences]
            Diff --> Final[ProcessResult]
        ```
        """
        remarks = []
        geom_thematic = make_valid(geom_thematic)
        if geom_preresult is None or geom_preresult.is_empty:
            geom_preresult = GeometryCollection()
            remark = ProcessRemark.RESULT_EMPTY_RETURNED
            remarks.append(remark)
            self.logger.feedback_warning(remark)

        elif to_multi(geom_preresult).geom_type != to_multi(geom_thematic).geom_type:
            geom_preresult = GeometryCollection()
            remark = ProcessRemark.CHANGED_GEOMETRYTYPE_EMPTY_RETURNED
            remarks.append(remark)
            self.logger.feedback_warning(remark)

        if geom_preresult.geom_type in [
            "Point",
            "MultiPoint",
            "LineString",
            "MultiLineString",
            "GeometryCollection",
        ]:
            result_diff_plus = safe_difference(
                geom_preresult,
                buffer_pos(geom_thematic, correction_distance),
            )
            result_diff_min = safe_difference(
                geom_thematic,
                buffer_pos(geom_preresult, correction_distance),
            )
            result_diff = safe_unary_union([result_diff_plus, result_diff_min])
            return union_process_result(
                {
                    "result": geom_preresult,
                    "result_diff": result_diff,
                    "result_diff_plus": result_diff_plus,
                    "result_diff_min": result_diff_min,
                    "result_relevant_intersection": relevant_intersection,
                    "result_relevant_diff": relevant_diff,
                    "properties": {REMARK_FIELD_NAME: remarks},
                }
            )

        buffer_distance = relevant_distance / 2
        result = []
        geom_thematic_for_add_delete = geom_thematic

        if self.config.od_strategy == OpenDomainStrategy.EXCLUDE:
            geom_thematic_for_add_delete = safe_intersection(
                geom_thematic_for_add_delete, reference_union
            )
            geom_preresult = safe_intersection(geom_preresult, reference_union)

        if not (geom_thematic is None or geom_thematic.is_empty):
            if (
                get_shape_index(geom_thematic.area, geom_thematic.length)
                > self.config.threshold_circle_ratio
            ):
                remark = ProcessRemark.INPUT_CIRCLE
                remarks.append(remark)
                self.logger.feedback_debug(remark)
                return union_process_result(
                    {
                        "result": geom_thematic,
                        "properties": {REMARK_FIELD_NAME: remarks},
                    }
                )

            if geometric_equality(
                geom_preresult,
                geom_thematic,
                correction_distance=correction_distance,
                mitre_limit=mitre_limit,
            ):
                remark = ProcessRemark.RESULT_UNCHANGED
                remarks.append(remark)
                self.logger.feedback_debug(remark)
                return union_process_result(
                    {
                        "result": geom_thematic,
                        "properties": {REMARK_FIELD_NAME: remarks},
                    }
                )

        geom_thematic_dissolved = buffer_pos(
            buffer_neg(
                buffer_pos(
                    geom_preresult, correction_distance, mitre_limit=mitre_limit
                ),
                2 * correction_distance,
                mitre_limit=mitre_limit,
            ),
            correction_distance,
            mitre_limit=mitre_limit,
        )

        geom_diff_add = safe_difference(
            geom_thematic_for_add_delete, geom_thematic_dissolved
        )
        geom_diff_delete = safe_difference(
            geom_thematic_dissolved, geom_thematic_for_add_delete
        )

        geom_diff_removed = safe_difference(
            geom_thematic_dissolved,
            safe_intersection(
                geom_diff_delete,
                buffer_neg_pos(
                    geom_diff_delete, buffer_distance, mitre_limit=mitre_limit
                ),
            ),
        )
        geom_diff_removed_added = safe_union(
            geom_diff_removed,
            safe_intersection(
                geom_diff_add,
                buffer_neg_pos(geom_diff_add, buffer_distance, mitre_limit=mitre_limit),
            ),
        )

        geom_thematic_preresult = buffer_pos(
            buffer_neg(
                buffer_pos(
                    geom_diff_removed_added,
                    correction_distance,
                    mitre_limit=mitre_limit,
                ),
                2 * correction_distance,
                mitre_limit=mitre_limit,
            ),
            correction_distance,
            mitre_limit=mitre_limit,
        )

        geom_thematic_cleaned_holes = fill_and_remove_gaps(
            geom_thematic_preresult, buffer_distance
        )

        geom_thematic_result = buffer_pos(
            buffer_neg(
                buffer_pos(
                    geom_thematic_cleaned_holes,
                    correction_distance,
                    mitre_limit=mitre_limit,
                ),
                2 * correction_distance,
                mitre_limit=mitre_limit,
            ),
            correction_distance,
            mitre_limit=mitre_limit,
        )
        geom_thematic_result = make_valid(remove_repeated_points(geom_thematic_result))

        if geom_thematic_result.is_empty or geom_thematic_result is None:
            geom_thematic_result = GeometryCollection()
            remark = ProcessRemark.RESULT_EMPTY_RETURNED
            remarks.append(remark)
            self.logger.feedback_warning(remark)

        result.append(geom_thematic_result)
        geom_thematic_result = safe_unary_union(result)

        geom_result_diff = buffer_pos(
            buffer_neg(
                safe_symmetric_difference(geom_thematic_result, geom_thematic),
                correction_distance,
                mitre_limit=mitre_limit,
            ),
            correction_distance,
            mitre_limit=mitre_limit,
        )
        geom_result_diff_plus = buffer_pos(
            buffer_neg(
                safe_difference(geom_thematic_result, geom_thematic),
                correction_distance,
                mitre_limit=mitre_limit,
            ),
            correction_distance,
            mitre_limit=mitre_limit,
        )
        geom_result_diff_min = buffer_pos(
            buffer_neg(
                safe_difference(geom_thematic, geom_thematic_result),
                correction_distance,
                mitre_limit=mitre_limit,
            ),
            correction_distance,
            mitre_limit=mitre_limit,
        )

        return union_process_result(
            {
                "result": geom_thematic_result,
                "result_diff": geom_result_diff,
                "result_diff_plus": geom_result_diff_plus,
                "result_diff_min": geom_result_diff_min,
                "result_relevant_intersection": relevant_intersection,
                "result_relevant_diff": relevant_diff,
                "properties": {REMARK_FIELD_NAME: remarks},
            }
        )

    def _create_virtual_reference(
        self,
        geometry: BaseGeometry,
        relevant_distance: float,
        reference_union: BaseGeometry,
        correction_distance: float,
        mitre_limit: float,
        use_outer_boundary: bool = False,
    ) -> BaseGeometry:
        """
        Creates a 'virtual reference polygon' for areas not covered by reference data.

        Useful for 'Onbekend Terrein' (Open Domain) strategies.

        Parameters
        ----------
        geometry : BaseGeometry
            The thematic geometry.
        relevant_distance : float
            Maximum alignment distance.
        reference_union : BaseGeometry
            Union of all actual reference features.
        correction_distance : float
            Small buffer for noise reduction.
        mitre_limit : float
            Mitre limit for buffering.
        use_outer_boundary : bool, optional
            If True, considers the outer boundary for the virtual polygon.

        Returns
        -------
        BaseGeometry
            The generated virtual reference geometry.
        """
        buffer_distance = relevant_distance / 2
        geom_thematic_buffered = make_valid(
            buffer_pos(
                geometry,
                (self.config.buffer_multiplication_factor * relevant_distance)
                + correction_distance,
                mitre_limit,
            )
        )
        clip_ref_thematic_buffered = safe_intersection(
            reference_union, geom_thematic_buffered
        )
        virtual_reference = safe_difference(
            geom_thematic_buffered, clip_ref_thematic_buffered
        )

        if use_outer_boundary:
            geom_1 = safe_difference(geometry, virtual_reference)
            geom_2 = buffer_neg_pos(geom_1, buffer_distance)
            geom_3 = safe_intersection(geom_2, geometry)
            virtual_reference = safe_unary_union([geom_3, virtual_reference])

        return virtual_reference

    @staticmethod
    def _add_multi_polygons_from_geom_to_array(
        geom: BaseGeometry, array: List[BaseGeometry]
    ) -> List[BaseGeometry]:
        """
        Extracts polygons from a geometry and appends them to an array.

        Parameters
        ----------
        geom : BaseGeometry
            Input geometry (could be a GeometryCollection).
        array : List[BaseGeometry]
            Target list to append polygons to.

        Returns
        -------
        List[BaseGeometry]
            The updated list containing extracted polygons.
        """
        if not (geom.is_empty or geom is None):
            geometry_collection = GeometryCollection(geom)
            for g in geometry_collection.geoms:
                g = make_valid(g)
                if str(g.geom_type) in ["Polygon", "MultiPolygon"]:
                    array.append(g)
        return array

    @staticmethod
    def _calculate_inner_outer(
        input_geometry: BaseGeometry, relevant_distance: float, max_outer_buffer: int
    ) -> tuple:
        """
        Splits a polygon into inner and outer zones for optimized processing.

        Parameters
        ----------
        input_geometry : BaseGeometry
            The geometry to split.
        relevant_distance : float
            The distance used to define the 'outer' shell.
        max_outer_buffer : int
            Value that is used to calculate the boundary of a thematic geometry wherefor the calculation has to be done. (inner part is added)

        Returns
        -------
        tuple
            A tuple of (inner_geometry, outer_geometry).
        """
        input_geometry = safe_unary_union(get_parts(input_geometry))
        input_geometry_inner = buffer_neg(input_geometry, relevant_distance)
        input_geometry_double_inner = buffer_neg(
            input_geometry, 2 * relevant_distance + max_outer_buffer
        )
        input_geometry_outer = safe_difference(
            input_geometry, input_geometry_double_inner
        )
        return input_geometry_inner, input_geometry_outer


class SnapGeometryProcessor(BaseProcessor):
    """
    Processor that aligns geometries by snapping them to the reference data.

    This processor uses the snapping algorithm to pull the vertices of the
    thematic geometry towards the nearest components of the reference data
    within a specified distance.

    Attributes
    ----------
    processor_id : ProcessorID
        The unique identifier for this processor (ProcessorID.SNAP).
    """

    processor_id = ProcessorID.SNAP

    def process(
        self,
        *,
        correction_distance: float,
        reference_data: AlignerFeatureCollection,
        input_geometry: BaseGeometry,
        mitre_limit: float,
        relevant_distance: float,
        **kwargs: Any,
    ) -> ProcessResult:
        """
        Aligns the input geometry by snapping its vertices to the reference geometries.

        The process considers the Open Domain (OD) strategy to determine how
        areas not covered by reference features should be handled (e.g.,
        ignored, kept as-is, or used as a virtual snapping target).

        Parameters
        ----------
        correction_distance : float
            Distance used for cleaning geometric noise.
        reference_data : AlignerFeatureCollection
            The collection of reference geometries.
        input_geometry : BaseGeometry
            The thematic geometry to be snapped.
        mitre_limit : float
            Mitre limit for buffering operations.
        relevant_distance : float
            The maximum distance within which snapping occurs.
        **kwargs : Any
            Additional arguments passed to the processor.

        Returns
        -------
        ProcessResult
            A dictionary containing the snapped result and difference metrics.

        Notes
        -----
        The snapping logic is influenced by the `od_strategy` (Open Domain):



        ```{mermaid}
        graph TD
            In[Input Geometry] --> OD{OD Strategy?}
            OD -- EXCLUDE --> Snap[Snap to Real Refs]
            OD -- AS_IS --> Keep[Keep OD part as-is]
            OD -- OTHER --> Virtual[Create Virtual Ref]
            Virtual --> SnapAll[Snap to Real + Virtual Refs]
            Keep --> Merge[Merge snapped & as-is parts]
            Snap --> Post[Post-process Result]
            SnapAll --> Post
            Merge --> Post
        ```
        """
        self.check_area_limit(input_geometry)
        snapped = []
        virtual_reference = Polygon()
        snap_strategy = self.config.snap_strategy

        # CALCULATE INNER and OUTER INPUT GEOMETRY for performance optimization on big geometries
        # combine all parts of the input geometry to one polygon
        input_geometry_inner, input_geometry_outer = self._calculate_inner_outer(
            input_geometry, relevant_distance, self.config.max_outer_buffer
        )
        # get a list of all ref_ids that are intersecting the thematic geometry; we take it bigger because we want to check if there are also reference geometries on a relevant distance.
        input_geometry_outer_buffered = buffer_pos(
            input_geometry_outer,
            relevant_distance * self.config.buffer_multiplication_factor,
        )
        ref_intersections = reference_data.items.take(
            reference_data.tree.query(input_geometry_outer_buffered)
        ).tolist()
        ref_intersections = reference_data.items.take(
            reference_data.tree.query(input_geometry_outer_buffered)
        ).tolist()

        ref_intersections_geoms = []
        for key_ref in ref_intersections:
            ref_geom = reference_data[key_ref].geometry
            ref_intersections_geoms.append(ref_geom)

        # Handle Open Domain (OD)logic
        if self.config.od_strategy != OpenDomainStrategy.EXCLUDE:
            virtual_reference = self._create_virtual_reference(
                input_geometry,
                relevant_distance,
                reference_data.union,
                correction_distance,
                mitre_limit,
                False,
            )

        if self.config.od_strategy == OpenDomainStrategy.EXCLUDE:
            pass
        elif self.config.od_strategy == OpenDomainStrategy.AS_IS:
            # Intersection with virtual reference is kept as original (no snapping)
            geom_od = safe_intersection(input_geometry, virtual_reference)
            snapped.append(geom_od)
        else:
            # Virtual reference is added to the list of snap targets
            ref_intersections_geoms.append(virtual_reference)

        # Execute the core snapping algorithm
        ref_geometrycollection = GeometryCollection(ref_intersections_geoms)
        snapped_geom = snap_geometry_to_reference(
            input_geometry,
            ref_geometrycollection,
            snap_strategy,
            self.config.snap_max_segment_length,
            relevant_distance,
        )
        snapped.append(snapped_geom)

        # Merge parts and clean the result
        geom_preresult = safe_unary_union(snapped)

        result_dict = self._postprocess_preresult(
            geom_preresult,
            input_geometry,
            GeometryCollection(),  # Relevant intersection is handled internally by snap
            GeometryCollection(),  # Relevant diff is handled internally by snap
            relevant_distance,
            reference_data.union,
            mitre_limit,
            correction_distance,
        )

        return result_dict


class DieussaertGeometryProcessor(BaseProcessor):
    """
    Processor implementing the Dieussaert area-based alignment algorithm.

    Unlike vertex snapping, this processor evaluates the overlap between thematic
    geometries and reference polygons. It decides per reference feature whether
    to include it fully, partially, or exclude it based on area-based thresholds.

    Attributes
    ----------
    processor_id : ProcessorID
        The unique identifier for this processor (ProcessorID.DIEUSSAERT).
    """

    processor_id = ProcessorID.DIEUSSAERT

    def process(
        self,
        *,
        input_geometry: BaseGeometry,
        reference_data: AlignerFeatureCollection,
        relevant_distance: float,
        mitre_limit: float,
        correction_distance: float,
        **kwargs: Any,
    ) -> ProcessResult:
        """
        Coordinates the alignment process for single or multi-geometries.

        If `multi_as_single_modus` is disabled, MultiPolygons are processed
        piecewise and merged afterward to ensure stability.

        Parameters
        ----------
        input_geometry : BaseGeometry
            The thematic geometry to align.
        reference_data : AlignerFeatureCollection
            The reference (target) dataset.
        relevant_distance : float
            The distance threshold for alignment decisions.
        mitre_limit : float
            Mitre limit for buffering operations.
        correction_distance : float
            Distance used for noise and sliver removal.
        **kwargs : Any
            Additional processor arguments.

        Returns
        -------
        ProcessResult
            The merged or single result of the alignment process.

        Notes
        -----
        The algorithm optimizes performance by splitting the input into an
        'inner' and 'outer' zone. Only the 'outer' zone (near the boundaries)
        is evaluated against reference data.



        ```{mermaid}
        graph TD
            Start[Input Geometry] --> Split[Split: Inner vs Outer]
            Split --> Query[Spatial Query Reference Data]
            Query --> OD[Process Open Domain]
            Query --> Ref[Process Intersecting Refs]
            OD --> Combine[Combine Results]
            Ref --> Combine
            Combine --> Post[Post-processing]
            Post --> End[ProcessResult]
        ```
        """
        if not isinstance(input_geometry, (Polygon, MultiPolygon)):
            raise ValueError(
                "Dieussaert algorithm can only be used when input geometry is polygon or multipolygon."
            )
        self.check_area_limit(input_geometry)
        if (
            not self.config.multi_as_single_modus
            or input_geometry is None
            or input_geometry.is_empty
            or len(to_multi(input_geometry).geoms) == 1
        ):
            return self._process(
                input_geometry=input_geometry,
                relevant_distance=relevant_distance,
                mitre_limit=mitre_limit,
                reference_data=reference_data,
                correction_distance=correction_distance,
            )
        else:
            input_geometry = to_multi(input_geometry)
            list_with_process_results = []
            for input_geom_single in input_geometry.geoms:
                process_result = self._process(
                    input_geometry=input_geom_single,
                    relevant_distance=relevant_distance,
                    mitre_limit=mitre_limit,
                    reference_data=reference_data,
                    correction_distance=correction_distance,
                )
                list_with_process_results.append(process_result)
            return self._merge_process_results(list_with_process_results)

    def _merge_process_results(
        self, list_process_results: list[ProcessResult]
    ) -> ProcessResult:
        """
        Merges multiple ProcessResult objects into a single result.

        This is used when a MultiPolygon is processed as individual components.
        It performs a unary union on all geometric keys and aggregates remarks.

        Parameters
        ----------
        list_process_results : list[ProcessResult]
            List of individual results to merge.

        Returns
        -------
        ProcessResult
            A single aggregated ProcessResult.
        """
        if len(list_process_results) == 1:
            return list_process_results[0]

        merged_process_result = ProcessResult()
        for process_result in list_process_results:
            for key in process_result:
                value = process_result[key]
                if key in ["metadata", "observation"]:
                    # This will add the needed keys to the processResult. The values are filled afterwards
                    if key in merged_process_result:
                        continue
                    else:
                        merged_process_result[key] = {}
                    continue
                if key == "properties":
                    if key in merged_process_result:
                        existing_remarks = merged_process_result[key][REMARK_FIELD_NAME]
                    else:
                        merged_process_result[key] = {}
                        existing_remarks = []

                    existing_remarks.extend(value[REMARK_FIELD_NAME])
                    merged_process_result[key][REMARK_FIELD_NAME] = existing_remarks
                    continue
                if isinstance(value, BaseGeometry):
                    geom = value
                    if geom is None:
                        geom = GeometryCollection()
                    if key in merged_process_result:
                        existing = merged_process_result[key]
                    else:
                        existing = GeometryCollection()
                    merged_process_result[key] = safe_unary_union([existing, geom])
                else:
                    raise ValueError(
                        f"Invalid element with key {str(key)} for ProcessResult"
                    )

        return merged_process_result

    def _process(
        self,
        *,
        input_geometry: BaseGeometry,
        reference_data: AlignerFeatureCollection,
        relevant_distance: float,
        mitre_limit: float,
        correction_distance: float,
    ) -> ProcessResult:
        """
        Internal core logic for the Dieussaert algorithm on a single geometry.
        """

        # CALCULATE INNER and OUTER INPUT GEOMETRY for performance optimization on big geometries
        # combine all parts of the input geometry to one polygon
        input_geometry_inner, input_geometry_outer = self._calculate_inner_outer(
            input_geometry, relevant_distance, self.config.max_outer_buffer
        )
        # get a list of all ref_ids that are intersecting the thematic geometry; we take it bigger because we want to check if there are also reference geometries on a relevant distance.
        input_geometry_outer_buffered = buffer_pos(
            input_geometry_outer,
            relevant_distance * self.config.buffer_multiplication_factor,
        )
        ref_intersections = reference_data.items.take(
            reference_data.tree.query(input_geometry_outer_buffered)
        ).tolist()

        ref_intersections_geoms = []
        for key_ref in ref_intersections:
            ref_geom = reference_data[key_ref].geometry
            ref_intersections_geoms.append(ref_geom)
            if not isinstance(ref_geom, (Polygon, MultiPolygon)):
                raise ValueError(
                    "Dieussaert algorithm can only be used when all reference geometries are polygons or multipolygons."
                )

        buffer_distance = relevant_distance / 2
        (
            preresult,
            relevant_intersection_array,
            relevant_diff_array,
        ) = self._calculate_intersection_between_geometry_and_od(
            input_geometry_outer,
            input_geometry_inner,
            relevant_distance,
            reference_data.union,
            mitre_limit,
            correction_distance,
        )

        for geom_reference in ref_intersections_geoms:
            geom_intersection = safe_intersection(input_geometry_outer, geom_reference)
            if geom_intersection.is_empty or geom_intersection is None:
                continue
            self.logger.feedback_debug("calculate intersection")
            (
                geom,
                relevant_intersection,
                relevant_diff,
            ) = self._calculate_geom_by_intersection_and_reference(
                geom_intersection=geom_intersection,
                geom_reference=geom_reference,
                input_geometry_inner=input_geometry_inner,
                is_open_domain=False,
                buffer_distance=buffer_distance,
                mitre_limit=mitre_limit,
            )
            self.logger.feedback_debug("intersection calculated")
            preresult = self._add_multi_polygons_from_geom_to_array(geom, preresult)
            relevant_intersection_array = self._add_multi_polygons_from_geom_to_array(
                relevant_intersection, relevant_intersection_array
            )
            relevant_diff_array = self._add_multi_polygons_from_geom_to_array(
                relevant_diff, relevant_diff_array
            )
        relevant_intersection = safe_unary_union(relevant_intersection_array)
        if relevant_intersection is None or relevant_intersection.is_empty:
            relevant_intersection = Polygon()
        relevant_diff = safe_unary_union(relevant_diff_array)
        if relevant_diff is None or relevant_diff.is_empty:
            relevant_diff = Polygon()
        preresult.append(input_geometry_inner)
        geom_preresult = safe_unary_union(preresult)
        process_result = self._postprocess_preresult(
            geom_preresult,
            input_geometry,
            relevant_intersection,
            relevant_diff,
            relevant_distance,
            reference_data.union,
            mitre_limit,
            correction_distance,
        )
        return process_result

    def _od_snap(
        self, geometry, relevant_distance, snap_strategy, mitre_limit, reference_union
    ):
        # all OD-parts wil be added AS IS
        geom_thematic_od = safe_difference(geometry, reference_union)
        if geom_thematic_od is None or geom_thematic_od.is_empty:
            return geom_thematic_od
        geom_thematic_od = to_multi(geom_thematic_od, geomtype="Polygon")
        out = []
        for p in geom_thematic_od.geoms:
            reference = safe_intersection(
                reference_union,
                make_valid(
                    buffer_pos(
                        p,
                        self.config.buffer_multiplication_factor * relevant_distance,
                        mitre_limit=mitre_limit,
                    )
                ),
            )
            p_snapped = snap_geometry_to_reference(
                p,
                reference,
                max_segment_length=self.config.partial_snap_max_segment_length,
                snap_strategy=snap_strategy,
                tolerance=relevant_distance,
            )
            out.append(p_snapped)
        return safe_unary_union(out)

    @staticmethod
    def _od_full_area(geometry, relevant_distance, reference_union, mitre_limit):
        buffer_distance = relevant_distance / 2
        geom_theme_od = safe_difference(geometry, reference_union)
        geom_theme_min_buffered = buffer_neg(
            buffer_pos(
                buffer_neg(
                    geometry,
                    relevant_distance,
                    mitre_limit=mitre_limit,
                ),
                buffer_distance,
                mitre_limit=mitre_limit,
            ),
            buffer_distance,
            mitre_limit=mitre_limit,
        )
        geom_theme_od_clipped_min_buffered = safe_intersection(
            geom_theme_min_buffered, geom_theme_od
        )
        geom_theme_od_min_clipped_plus_buffered = buffer_pos(
            geom_theme_od_clipped_min_buffered,
            relevant_distance,
            mitre_limit=mitre_limit,
        )
        geom_theme_od_min_clipped_plus_buffered_clipped = safe_intersection(
            geom_theme_od_min_clipped_plus_buffered, geom_theme_od
        )
        geom_thematic_od = geom_theme_od_min_clipped_plus_buffered_clipped
        return geom_thematic_od

    def _od_snap_all_side(
        self,
        geometry,
        input_geometry_inner,
        relevant_distance,
        correction_distance,
        mitre_limit,
        reference_union,
        outer=False,
    ):
        """

        :param geometry:
        :param input_geometry_inner:
        :param relevant_distance:
        :param outer: when outer is True, the outer boundary is used, inner is not used
        :return:
        """
        buffer_distance = relevant_distance / 2
        relevant_difference_array = []
        relevant_intersection_array = []
        geom_reference = self._create_virtual_reference(
            geometry,
            relevant_distance,
            reference_union,
            correction_distance,
            mitre_limit,
            outer,
        )

        geom_thematic = geometry
        if geom_reference.is_empty or geom_reference is None:
            geom_thematic_od = geom_reference
        else:
            geom_intersection = safe_intersection(geom_reference, geom_thematic)
            # geom_intersection = safe_unary_union([geom_intersection])

            # calculate the geometry for the 'virtual' OD parcel
            (
                geom_thematic_od,
                geom_relevant_intersection,
                geom_relevant_diff,
            ) = self._calculate_geom_by_intersection_and_reference(
                geom_intersection=geom_intersection,
                geom_reference=geom_reference,
                input_geometry_inner=input_geometry_inner,
                is_open_domain=True,
                buffer_distance=buffer_distance,
                mitre_limit=mitre_limit,
            )

            relevant_intersection_array = self._add_multi_polygons_from_geom_to_array(
                geom_relevant_intersection, []
            )
            relevant_difference_array = self._add_multi_polygons_from_geom_to_array(
                geom_relevant_diff, []
            )
            # geom_thematic_od = safe_unary_union(geom_thematic_od_array)
        return geom_thematic_od, relevant_difference_array, relevant_intersection_array

    def _calculate_intersection_between_geometry_and_od(
        self,
        input_geometry,
        input_geometry_inner,
        relevant_distance,
        reference_union,
        mitre_limit,
        correction_distance,
    ):
        """
        Calculates the intersecting parts between a thematic geometry and the Open Domain( domain, not covered by reference-polygons)
        :param input_geometry:
        :param relevant_distance:
        :return:
        """
        # Calculate the intersection between thematic and Open Domain
        # buffer_distance = relevant_distance / 2
        relevant_intersection_array = []
        relevant_difference_array = []
        geom_thematic_od = Polygon()

        if self.config.od_strategy == OpenDomainStrategy.AS_IS:
            # All parts that are not covered by the reference layer are added to the
            #         resulting geometry AS IS
            self.logger.feedback_debug("OD-strategy AS IS")
            # all OD-parts wil be added AS IS
            geom_thematic_od = safe_difference(input_geometry, reference_union)

        elif (
            self.config.od_strategy == OpenDomainStrategy.SNAP_INNER_SIDE
            or self.config.od_strategy == OpenDomainStrategy.EXCLUDE
        ):
            # integrates the entire inner area of the input geometry,
            # so Open Domain of the inner area is included in the result
            self.logger.feedback_debug("OD-strategy OD_SNAP_INNER_SIDE or EXCLUDE")
            geom_thematic_od = self._od_full_area(
                input_geometry, relevant_distance, reference_union, mitre_limit
            )

        elif self.config.od_strategy == OpenDomainStrategy.SNAP_ALL_SIDE:
            #  Inner & Outer-reference-boundaries are used.
            # integrates the entire inner area of the input geometry,
            self.logger.feedback_debug("OD-strategy OD-SNAP_ALL_SIDE")
            (
                geom_thematic_od,
                relevant_difference_array,
                relevant_intersection_array,
            ) = self._od_snap_all_side(
                input_geometry,
                input_geometry_inner,
                relevant_distance,
                correction_distance,
                mitre_limit,
                reference_union,
                outer=False,
            )
            # This part calculates the full area
            geom_theme_od_min_clipped_plus_buffered_clipped = self._od_full_area(
                geometry=input_geometry,
                relevant_distance=relevant_distance,
                reference_union=reference_union,
                mitre_limit=mitre_limit,
            )
            # UNION of both elements
            geom_thematic_od = safe_union(
                geom_theme_od_min_clipped_plus_buffered_clipped, geom_thematic_od
            )
        elif self.config.od_strategy == OpenDomainStrategy.SNAP_PREFER_VERTICES:
            self.logger.feedback_debug("OD-strategy SNAP_PREFER_VERTICES")
            geom_thematic_od = self._od_snap(
                geometry=input_geometry,
                relevant_distance=relevant_distance,
                snap_strategy=SnapStrategy.PREFER_VERTICES,
                mitre_limit=mitre_limit,
                reference_union=reference_union,
            )

        elif self.config.od_strategy == OpenDomainStrategy.SNAP_NO_PREFERENCE:
            self.logger.feedback_debug("OD-strategy SNAP_NO_PREFERENCE")
            geom_thematic_od = self._od_snap(
                geometry=input_geometry,
                relevant_distance=relevant_distance,
                snap_strategy=SnapStrategy.NO_PREFERENCE,
                mitre_limit=mitre_limit,
                reference_union=reference_union,
            )

        elif self.config.od_strategy == OpenDomainStrategy.SNAP_ONLY_VERTICES:
            self.logger.feedback_debug("OD-strategy SNAP_ONLY_VERTICES")
            geom_thematic_od = self._od_snap(
                geometry=input_geometry,
                relevant_distance=relevant_distance,
                snap_strategy=SnapStrategy.ONLY_VERTICES,
                mitre_limit=mitre_limit,
                reference_union=reference_union,
            )

        # ADD THEMATIC_OD
        preresult = self._add_multi_polygons_from_geom_to_array(geom_thematic_od, [])
        return (
            preresult,
            relevant_intersection_array,
            relevant_difference_array,
        )

    def _calculate_geom_by_intersection_and_reference(
        self,
        *,
        geom_intersection: BaseGeometry,
        geom_reference: BaseGeometry,
        input_geometry_inner: BaseGeometry,
        is_open_domain: bool,
        buffer_distance: float,
        mitre_limit: float,
    ) -> tuple[BaseGeometry, BaseGeometry, BaseGeometry]:
        """
        Decides the resulting geometry for a specific intersection area.

        This is the heart of the area-based decision logic. It calculates the
        percentage of overlap and applies thresholds for inclusion/exclusion.

        Parameters
        ----------
        geom_intersection : BaseGeometry
            The intersection between input and reference.
        geom_reference : BaseGeometry
            The specific reference feature being evaluated.
        input_geometry_inner : BaseGeometry
            The protected inner part of the input.
        is_open_domain : bool
            Whether this calculation is for an area without reference features.
        buffer_distance : float
            Buffer used to define "relevance".
        mitre_limit : float
            Mitre limit for buffering.

        Returns
        -------
        tuple
            A tuple containing (result_geometry, relevant_intersection, relevant_difference).

        Notes
        -----
        Overlap criteria:
        -   **Full Inclusion**: Overlap > `threshold_inclusion_percentage`.
        -   **Full Exclusion**: Overlap < `threshold_exclusion_percentage`.
        -   **Partial Alignment**: Calculated via negative/positive buffering
            when overlap falls between thresholds.


        """
        od_overlap = 111  # define a specific value for defining overlap of OD
        if geom_reference.area == 0:
            overlap = od_overlap  # Open Domain

        else:
            overlap = geom_intersection.area * 100 / geom_reference.area

        if (
            overlap < self.config.threshold_exclusion_percentage
            or geom_intersection.area < self.config.threshold_exclusion_area
        ):
            return Polygon(), Polygon(), Polygon()

        if (
            overlap >= self.config.threshold_inclusion_percentage
            and not overlap == od_overlap
        ):
            return geom_reference, geom_reference, Polygon()

        geom_difference = safe_difference(geom_reference, geom_intersection)
        geom_relevant_intersection = buffer_neg(
            geom_intersection,
            buffer_distance,
            mitre_limit=mitre_limit,
        )
        geom_intersection_inner = safe_intersection(
            geom_intersection, input_geometry_inner
        )  # this part is the intersection-part that always has to be kept, because it is inside the inner_geometry

        geom_relevant_difference = buffer_neg(
            geom_difference,
            buffer_distance,
            mitre_limit=mitre_limit,
        )

        if (
            not geom_intersection_inner.is_empty
            and geom_relevant_intersection.is_empty
            and not geom_relevant_difference.is_empty
        ):
            geom_x = safe_intersection(
                geom_intersection,
                buffer_pos(geom_intersection_inner, 2 * buffer_distance),
            )

            geom_x = snap_geometry_to_reference(
                geom_x,
                geom_reference,
                max_segment_length=self.config.partial_snap_max_segment_length,
                snap_strategy=self.config.partial_snap_strategy,
                tolerance=2 * buffer_distance,
            )

            geom = geom_x
        elif (
            not geom_relevant_intersection.is_empty
            and not geom_relevant_difference.is_empty
        ):
            # relevant intersection and relevant difference

            geom_x = safe_difference(
                geom_reference,
                safe_intersection(
                    geom_difference,
                    buffer_neg_pos(
                        geom_difference,
                        buffer_distance,
                        mitre_limit=mitre_limit,
                    ),
                ),
            )
            geom_x = buffer_neg_pos(geom_x, buffer_distance, mitre_limit=mitre_limit)

            geom_intersection_buffered = buffer_pos(
                geom_intersection, 2 * buffer_distance
            )
            geom_difference_2 = safe_difference(
                geom_reference, geom_intersection_buffered
            )
            geom_difference_2_buffered = buffer_pos(
                geom_difference_2, 2 * buffer_distance
            )

            geom_x = safe_difference(geom_x, geom_difference_2_buffered)

            geom_x = safe_intersection(geom_x, geom_reference)

            if self.config.partial_snapping:
                geom_x = snap_geometry_to_reference(
                    geom_x,
                    geom_reference,
                    max_segment_length=self.config.partial_snap_max_segment_length,
                    snap_strategy=self.config.partial_snap_strategy,
                    tolerance=2 * buffer_distance,
                )
            geom = safe_unary_union(
                [geom_x, geom_relevant_intersection, geom_intersection_inner]
            )

            # when calculating for OD, we create a 'virtual parcel'. When calculating this
            # virtual parcel, it is buffered to take outer boundaries into account.
            # This results in a side effect that there are extra non-logical parts included
            # in the result. The function below tries to exclude these non-logical parts.
            # see eo_id 206363 with relevant distance=0.2m and SNAP_ALL_SIDE
            if is_open_domain:
                geom = get_relevant_polygons_from_geom(
                    geom, buffer_distance, mitre_limit
                )
        elif (
            not geom_relevant_intersection.is_empty
            and geom_relevant_difference.is_empty
        ):
            geom = geom_reference
        elif (
            geom_relevant_intersection.is_empty
            and not geom_relevant_difference.is_empty
        ):
            geom = geom_relevant_intersection  # (=empty geometry)
        else:
            # No relevant intersection and no relevant difference
            if is_open_domain:
                # geom = geom_relevant_intersection  # (=empty geometry)
                geom = snap_geometry_to_reference(
                    geom_intersection,
                    geom_reference,
                    snap_strategy=self.config.partial_snap_strategy,
                    tolerance=2 * buffer_distance,
                    max_segment_length=self.config.partial_snap_max_segment_length,
                )
            elif not geom_intersection_inner.is_empty:
                geom_intersection_buffered = buffer_pos(
                    geom_intersection, 2 * buffer_distance
                )
                geom_difference_2 = safe_difference(
                    geom_reference, geom_intersection_buffered
                )
                geom_difference_2_buffered = buffer_pos(
                    geom_difference_2, 2 * buffer_distance
                )
                geom = safe_difference(geom_reference, geom_difference_2_buffered)
            elif self.config.threshold_overlap_percentage < 0:
                # if we take a value of -1, the original border will be used
                geom = geom_intersection
            elif overlap > self.config.threshold_overlap_percentage:
                geom = geom_reference
            else:
                geom = geom_relevant_intersection  # (=empty geometry)
        return geom, geom_relevant_intersection, geom_relevant_difference


class NetworkGeometryProcessor(BaseProcessor):
    """
    Processor that aligns geometries based on a linear network.

    This processor decomposes polygons into their exterior and interior linear
    rings and aligns these boundaries to the linear elements (lines and points)
    found in the reference dataset.

    Attributes
    ----------
    processor_id : ProcessorID
        The unique identifier for this processor (ProcessorID.NETWORK).
    """

    processor_id = ProcessorID.NETWORK

    def process(
        self,
        *,
        input_geometry: BaseGeometry,
        reference_data: AlignerFeatureCollection,
        mitre_limit: float,
        correction_distance: float,
        relevant_distance: float,
        **kwargs: Any,
    ) -> ProcessResult:
        """
        Aligns the boundaries of the input geometry to a reference network.

        The process buffers the input to find relevant network elements,
        processes exterior and interior rings separately, and reconstructs
        the polygon after alignment.

        Parameters
        ----------
        input_geometry : BaseGeometry
            The thematic geometry (Polygon or MultiPolygon) to align.
        reference_data : AlignerFeatureCollection
            The reference dataset, specifically using its `elements` property.
        mitre_limit : float
            Mitre limit for buffering operations.
        correction_distance : float
            Distance used for cleaning and noise reduction.
        relevant_distance : float
            The maximum distance to search for network elements.
        **kwargs : Any
            Additional arguments passed to the processor.

        Returns
        -------
        ProcessResult
            A dictionary containing the reconstructed and cleaned polygon.

        Notes
        -----
        The network processing follows a "deconstruct-align-reconstruct" flow:



        ```{mermaid}
        graph TD
            In[Input Polygon] --> Decon[Deconstruct: Exterior & Interiors]
            Decon --> Buff[Buffer Input to Find Network]
            Buff --> Align[Align Segments to Network Elements]
            Align --> Recon[Reconstruct Polygon Rings]
            Recon --> Post[Post-processing & Sliver Removal]
            Post --> End[Final ProcessResult]
        ```
        """
        self.check_area_limit(input_geometry)
        input_geometry = to_multi(input_geometry)

        # Determine the search area for relevant network elements
        input_geometry_buffered = buffer_pos(
            input_geometry,
            relevant_distance * self.config.buffer_multiplication_factor,
        )

        # Fetch linear/point elements from reference that fall within the buffer
        reference = safe_unary_union(
            safe_intersection(reference_data.elements, input_geometry_buffered)
        )

        geom_processed_list = []

        if isinstance(input_geometry, MultiPolygon):
            # Cast to MultiPolygon for consistent iteration
            for polygon in input_geometry.geoms:
                # 1. Process the outer boundary
                exterior = polygon.exterior
                exterior_processed = self._process_by_network(
                    exterior,
                    reference,
                    relevant_distance,
                    close_output=True,
                )

                # 2. Process all inner holes
                interiors = polygon.interiors
                interiors_processed = []
                for i in interiors:
                    i_processed = self._process_by_network(
                        i,
                        reference,
                        relevant_distance,
                        close_output=True,
                    )
                    interiors_processed.append(i_processed)

                # 3. Reconstruct the polygon
                geom_processed = Polygon(exterior_processed, interiors_processed)
                geom_processed_list.append(geom_processed)
        else:
            # Handling for non-polygonal geometries (e.g. LineStrings)
            for geom in input_geometry.geoms:
                geom_processed = self._process_by_network(
                    geom,
                    reference,
                    relevant_distance,
                    close_output=False,
                )
                geom_processed_list.append(geom_processed)

        # Merge all processed parts
        geom_processed = safe_unary_union(geom_processed_list)

        # Standard cleaning pipeline
        return self._postprocess_preresult(
            geom_processed,
            input_geometry,
            GeometryCollection(),
            GeometryCollection(),
            relevant_distance,
            reference_data.union,
            mitre_limit,
            correction_distance,
        )

    def _process_by_network(
        self,
        geom_to_process,
        reference,
        relevant_distance,
        close_output=False,
    ):
        buffer_distance = relevant_distance / 2
        geom_to_process_buffered = buffer_pos(geom_to_process, buffer_distance)
        reference_buffered = buffer_pos(reference, buffer_distance)
        overlap = safe_intersection(geom_to_process_buffered, reference_buffered)
        if overlap.is_empty:
            return geom_to_process
        overlap_buffered = buffer_pos(overlap, buffer_distance)
        reference_intersection = safe_intersection(reference, overlap_buffered)
        reference_intersection = safe_unary_union(reference_intersection)
        reference_intersection = to_multi(reference_intersection)
        if reference_intersection.is_empty:
            return geom_to_process
        thematic_difference = safe_difference(geom_to_process, overlap_buffered)
        thematic_difference = safe_unary_union(thematic_difference)
        thematic_difference = to_multi(thematic_difference)
        if self.config.snap_strategy != SnapStrategy.NO_PREFERENCE:
            reference_coords = MultiPoint(list(get_coords_from_geometry(reference)))
            reference_coords_intersection = to_multi(
                safe_intersection(reference_coords, overlap_buffered)
            )
        else:
            reference_coords_intersection = None
        if isinstance(geom_to_process, Point):
            p1, p2 = nearest_points(geom_to_process, reference_intersection)
            if (
                not reference_coords_intersection is None
                and not reference_coords_intersection.is_empty
            ):
                p1_vertices, p2_vertices = nearest_points(
                    geom_to_process, reference_coords_intersection
                )
                if not p2_vertices is None and not p2_vertices.is_empty:
                    p2 = p2_vertices
            geom_processed = p2
        else:
            geom_processed = self._get_processed_network_path(
                input_geometry=geom_to_process,
                reference_intersection=reference_intersection,
                thematic_difference=thematic_difference,
                relevant_distance=relevant_distance,
            )
            if (
                close_output
                and geom_processed is not None
                and not geom_processed.is_ring
            ):
                closed_coords = list(geom_processed.coords) + [geom_processed.coords[0]]
                geom_processed = LineString(closed_coords)
        return make_valid(geom_processed)

    def _get_thematic_points(self, geom_to_process, reference_intersection):
        """
        returns a MultiPoint geometry with points on the thematic geometry that are used for consistency while creating connection_lines
        :param geom_to_process:
        :param reference_intersection:
        :return: MultiPoint
        """
        geom_to_process_line = geom_to_process
        if isinstance(geom_to_process_line, LinearRing):
            geom_to_process_line = LineString(geom_to_process_line.coords)

        # geom_to_process_segmentized = segmentize(
        #     geom_to_process_line, self.config.partial_snap_max_segment_length
        # )
        geom_to_process_segmentized = geom_to_process_line
        # Split the line at all intersection points with the MultiLineString
        splitter = safe_unary_union(reference_intersection)
        try:
            geom_to_process_splitted = split(geom_to_process_segmentized, splitter)
        except (GeometryTypeError, ValueError):
            geom_to_process_splitted = geom_to_process_segmentized
        thematic_points = MultiPoint(
            list(get_coords_from_geometry(geom_to_process_splitted))
        )
        return thematic_points

    def _get_processed_network_path(
        self,
        input_geometry,
        reference_intersection,
        thematic_difference,
        relevant_distance,
    ):
        ref_multilinestring, ref_points = (
            self._multilinestring_multipoint_from_reference_intersection(
                reference_intersection
            )
        )
        thematic_points = self._get_thematic_points(
            input_geometry, reference_intersection
        )
        graph = build_custom_network(
            theme_multiline=thematic_difference,
            ref_multiline=ref_multilinestring,
            ref_points=ref_points,
            theme_points=thematic_points,
            gap_threshold=0.1,
            relevant_distance=relevant_distance,
        )
        return find_best_path_in_network(
            input_geometry, graph, self.config.snap_strategy, relevant_distance
        )

    def _multilinestring_multipoint_from_reference_intersection(
        self, reference_intersection
    ):
        if isinstance(reference_intersection, (LineString, MultiLineString)):
            return to_multi(reference_intersection), MultiPoint()
        if isinstance(reference_intersection, (Point, MultiPoint)):
            return MultiLineString(), to_multi(reference_intersection)
        if isinstance(reference_intersection, GeometryCollection):
            points = []
            lines = []
            for geom in reference_intersection.geoms:

                if isinstance(geom, (Point, MultiPoint)):
                    points.append(geom)
                elif isinstance(geom, (LineString, MultiLineString)):
                    lines.append(geom)
                elif isinstance(geom, (Polygon, MultiPolygon)):
                    lines.append(geom.boundary)
                else:
                    TypeError("Geometrytype not valid at this stage")
            return to_multi(safe_unary_union(lines)), to_multi(safe_unary_union(points))

        if isinstance(reference_intersection, (Polygon, MultiPolygon)):
            return to_multi(reference_intersection.boundary), MultiPoint()

        raise TypeError(
            "Reference could not be interpreted by NetworkGeometryProcessor"
        )


class AlignerGeometryProcessor(BaseProcessor):
    """
    Processor responsible for aligning thematic geometries to reference data.

    This class identifies the geometry type and delegates the alignment logic
    to specialized processors (Dieussaert or Network-based) while handling
    validation and edge cases like area limits and zero-distance processing.

    Attributes
    ----------
    processor_id : ProcessorID
        Unique identifier for the aligner processor.
    """

    processor_id = ProcessorID.ALIGNER

    def process(
        self,
        *,
        correction_distance,
        reference_data,
        input_geometry: BaseGeometry,
        mitre_limit,
        relevant_distance: float = 1.0,
        **kwargs,
    ) -> ProcessResult:
        """
        Process and align a single geometry based on reference data.

        The method validates the input geometry, checks against area constraints,
        and selects the appropriate sub-processor (Dieussaert for polygons or
        Network for linear/complex structures).

        Parameters
        ----------
        correction_distance : float
            The maximum distance a vertex can be moved during alignment.
        reference_data : GeoDataFrame or list[BaseGeometry]
            The target geometries to which the input_geometry should align.
        input_geometry : BaseGeometry
            The geometry to be processed. Supports Polygon, MultiPolygon,
            and linear geometries.
        mitre_limit : float
            The limit used for miter joins to prevent sharp spikes in corners.
        relevant_distance : float, default 1.0
            The search radius used to find nearby reference geometries.
            If set to 0, no processing is performed.
        **kwargs : dict
            Additional keyword arguments passed to sub-processors.

        Returns
        -------
        ProcessResult
            An object containing the aligned geometry and processing metadata
            (e.g., remarks or error logs).

        Raises
        ------
        ValueError
            If the input_geometry is a GeometryCollection.
            If the input_geometry exceeds the configured area_limit.

        Notes
        -----
        The logic follows this decision flow:

        ```{mermaid}
        graph TD
              A[Start Process] --> B{GeometryCollection?}
              B -- Yes --> C[Raise ValueError]
              B -- No --> D{Polygon/MultiPolygon?}
              D -- No --> E[Network Processor]
              D -- Yes --> F{Area > Limit?}
              F -- Yes --> G[Raise ValueError]
              F -- No --> H{RD == 0?}
              H -- Yes --> I[Return Original]
              H -- No --> J[Dieussaert Processor]
              J -- Success --> K[Return Result]
              J -- Fail --> E
              E --> K
        ```
        """

        if isinstance(input_geometry, GeometryCollection):
            raise ValueError(
                "GeometryCollection as input is not supported. Please use the individual geometries from the GeometryCollection as input."
            )
        elif isinstance(input_geometry, (Polygon, MultiPolygon)):
            # Processing thematic polygons
            self.check_area_limit(input_geometry)
            self.logger.feedback_debug("process geometry")

            # For calculations with RD=0 the original input is returned
            if relevant_distance == 0:
                remark = ProcessRemark.RESULT_UNCHANGED
                self.logger.feedback_debug(remark)
                remarks = [remark]
                process_result = ProcessResult()
                process_result["result"] = input_geometry
                process_result["properties"] = {REMARK_FIELD_NAME: remarks}
                return union_process_result(process_result)

            try:
                processor = DieussaertGeometryProcessor(
                    self.config,
                    self.logger.feedback,
                )
                return processor.process(
                    input_geometry=input_geometry,
                    reference_data=reference_data,
                    relevant_distance=relevant_distance,
                    mitre_limit=mitre_limit,
                    correction_distance=correction_distance,
                )
            except ValueError as e:
                self.logger.feedback_debug(
                    f"Dieussaert processing failed with error: {str(e)}. Trying Network-based processing."
                )
        processor = NetworkGeometryProcessor(self.config, self.logger.feedback)
        return processor.process(
            input_geometry=input_geometry,
            reference_data=reference_data,
            relevant_distance=relevant_distance,
            mitre_limit=mitre_limit,
            correction_distance=correction_distance,
        )


class TopologyProcessor(BaseProcessor):
    """
    Processor that aligns geometries while preserving topological relationships.

    Instead of processing features independently, the TopologyProcessor
    decomposes thematic data into unique 'arcs' (shared boundaries). Each arc
    is aligned once using network-based processing, ensuring that shared
    boundaries remain perfectly snapped together in the final output.

    Attributes
    ----------
    processor_id : ProcessorID
        The unique identifier for this processor (ProcessorID.TOPOLOGY).
    thematic_data : AlignerFeatureCollection, optional
        Cached reference to the full thematic dataset used for topology building.
    topo_thematic : dict, optional
        The generated TopoJSON-like structure of the thematic data.
    thematic_geometries_to_process : dict, optional
        Mapping of arc IDs to their respective LineString geometries.
    id_to_arcs : dict, optional
        Mapping of feature IDs to the list of arc IDs that form their boundary.
    wkb_to_id : dict, optional
        Reverse index mapping geometry WKB strings to feature IDs.
    """

    processor_id = ProcessorID.TOPOLOGY

    def __init__(self, config: ProcessorConfig, feedback: Any = None):
        """
        Initializes the TopologyProcessor with internal cache storage.
        """
        super().__init__(config, feedback)
        self.thematic_data = None
        self.topo_thematic = None
        self.thematic_geometries_to_process = None
        self.id_to_arcs = None
        self.wkb_to_id = None

    def process(
        self,
        *,
        input_geometry: BaseGeometry,
        reference_data: AlignerFeatureCollection,
        mitre_limit: float,
        correction_distance: float,
        relevant_distance: float,
        thematic_data: AlignerFeatureCollection,
        **kwargs: Any,
    ) -> ProcessResult:
        """
        Aligns a geometry by processing its topological arcs.

        This method identifies which unique arcs belong to the input geometry,
        aligns those arcs using the `NetworkGeometryProcessor`, and then
        reconstructs the final geometry.

        Parameters
        ----------
        input_geometry : BaseGeometry
            The specific thematic geometry to align.
        reference_data : AlignerFeatureCollection
            The reference dataset used for alignment.
        mitre_limit : float
            Mitre limit for buffering operations.
        correction_distance : float
            Distance used for cleaning and noise reduction.
        relevant_distance : float
            The maximum distance for alignment.
        thematic_data : AlignerFeatureCollection
            The full thematic collection required to build/query the topology.
        **kwargs : Any
            Additional arguments.

        Returns
        -------
        ProcessResult
            The reconstructed geometry where shared boundaries are consistently aligned.

        Notes
        -----
        The topological workflow ensures "gapless" alignment:



        ```{mermaid}
        graph TD
            In[Input Geometry] --> Cache{Cache Built?}
            Cache -- No --> Build[Build Topology: Extract Arcs]
            Build --> Cache
            Cache -- Yes --> Map[Map Geometry to Arc IDs]
            Map --> Loop[For each unique Arc]
            Loop --> Net[Align Arc via NetworkProcessor]
            Net --> Loop
            Loop --> Dissolve[Reconstruct Polygon from Aligned Arcs]
            Dissolve --> End[Final ProcessResult]
        ```
        """
        self.check_area_limit(input_geometry)
        self._build_topo_cache(thematic_data)

        # Identify the feature ID via its WKB
        id_thematic = self.wkb_to_id[input_geometry.wkb]

        # Get all unique arcs that make up this feature
        arcs_to_process = flatten_iter(self.id_to_arcs[id_thematic])
        processor = NetworkGeometryProcessor(config=self.config, feedback=self.feedback)
        process_results = {}

        # Align each arc individually (shared arcs are processed once per batch)
        for key in arcs_to_process:
            key = abs(key)
            geometry = self.thematic_geometries_to_process[key]
            process_results[key] = {}
            process_results[key][relevant_distance] = processor.process(
                correction_distance=correction_distance,
                mitre_limit=mitre_limit,
                reference_data=reference_data,
                input_geometry=geometry,
                relevant_distance=relevant_distance,
            )

        # Reconstruct the original polygon structure using the new arc geometries
        return _dissolve_topo(
            id_thematic,
            process_results,
            input_geometry,
            self.thematic_geometries_to_process,
            self.topo_thematic,
            relevant_distance,
        )

    def _build_topo_cache(self, thematic_data: AlignerFeatureCollection):
        """
        Builds and caches the topological structure of the thematic data.

        This is an expensive operation performed only once per dataset change.
        It decomposes all features into a set of non-overlapping arcs.

        Parameters
        ----------
        thematic_data : AlignerFeatureCollection
            The dataset to analyze for shared boundaries.

        Notes
        -----
        The cache consists of:
        1. A TopoJSON-like structure (`topo_thematic`).
        2. A mapping of arc IDs to LineStrings (`thematic_geometries_to_process`).
        3. A reverse index for quick feature lookup.
        """
        if self.thematic_data is None or thematic_data != self.thematic_data:
            self.thematic_geometries_to_process, self.topo_thematic = _generate_topo(
                thematic_data
            )
            self.wkb_to_id = build_reverse_index_wkb(
                {key: feat.geometry for key, feat in thematic_data.features.items()}
            )
            self.id_to_arcs = _topojson_id_to_arcs(self.topo_thematic)
            self.thematic_data = thematic_data

        return
