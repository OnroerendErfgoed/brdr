import heapq
import json
import os
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import wait
from datetime import datetime
from typing import Iterable

import numpy as np
from shapely import (
    GeometryCollection,
    Point,
    MultiPolygon,
    MultiPoint,
    segmentize,
    shortest_line,
    LinearRing,
)
from shapely import Polygon
from shapely import STRtree
from shapely import get_parts
from shapely import make_valid
from shapely import remove_repeated_points
from shapely import to_geojson
from shapely.errors import GeometryTypeError
from shapely.geometry.base import BaseGeometry
from shapely.geometry.linestring import LineString
from shapely.ops import nearest_points, split

from brdr import __version__
from brdr.constants import (
    DEFAULT_CRS,
    PREDICTION_SCORE,
    PREDICTION_COUNT,
    MAX_OUTER_BUFFER,
    PARTIAL_SNAP_MAX_SEGMENT_LENGTH,
    PARTIAL_SNAP_STRATEGY,
    PARTIAL_SNAPPING,
    RELEVANT_DISTANCE_DECIMALS,
    ID_THEME_FIELD_NAME,
    ID_REFERENCE_FIELD_NAME,
    SNAP_STRATEGY,
    SNAP_MAX_SEGMENT_LENGTH,
    BUFFER_MULTIPLICATION_FACTOR,
    DIFF_METRIC,
)
from brdr.constants import (
    LAST_VERSION_DATE,
    VERSION_DATE,
    DATE_FORMAT,
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
    AlignerResultType,
    AlignerInputType,
    SnapStrategy,
    FullStrategy,
    DiffMetric,
)
from brdr.geometry_utils import (
    buffer_neg,
    safe_unary_union,
    get_shape_index,
    geometric_equality,
    to_multi,
    snap_geometry_to_reference,
    extract_points_lines_from_geometry,
    shortest_connections_between_geometries,
    get_coords_from_geometry,
    find_best_path_in_network,
    prepare_network,
)
from brdr.geometry_utils import buffer_neg_pos
from brdr.geometry_utils import buffer_pos
from brdr.geometry_utils import fill_and_remove_gaps
from brdr.geometry_utils import safe_difference
from brdr.geometry_utils import safe_intersection
from brdr.geometry_utils import safe_symmetric_difference
from brdr.geometry_utils import safe_union
from brdr.loader import Loader
from brdr.logger import Logger
from brdr.topo import dissolve_topo, generate_topo
from brdr.typings import ProcessResult
from brdr.utils import (
    diffs_from_dict_processresults,
    multi_to_singles,
    is_brdr_formula,
)
from brdr.utils import geojson_from_dict
from brdr.utils import get_breakpoints_zerostreak
from brdr.utils import get_dict_geojsons_from_series_dict
from brdr.utils import merge_process_results
from brdr.utils import write_geojson


###################


class Aligner:
    """
    This class is used to compare and align the thematic data with the reference data.
    The reference data can be loaded in different ways, for example by using the GRB
    data.
    The thematic data can be loaded by using different Loaders: DictLoader, GeojsonLoader,...
    The class can be used to compare and aligne the thematic data with the reference data.
    """

    def __init__(
        self,
        *,
        feedback=None,
        relevant_distance=1,
        relevant_distances=[
            round(k, RELEVANT_DISTANCE_DECIMALS)
            for k in np.arange(0, 310, 10, dtype=int) / 100
        ],
        threshold_overlap_percentage=50,
        od_strategy=OpenDomainStrategy.SNAP_ALL_SIDE,
        crs=DEFAULT_CRS,
        multi_as_single_modus=True,
        preserve_topology=False,
        snap_strategy=SNAP_STRATEGY,
        snap_max_segment_length=SNAP_MAX_SEGMENT_LENGTH,
        partial_snapping=PARTIAL_SNAPPING,
        partial_snap_strategy=PARTIAL_SNAP_STRATEGY,
        partial_snap_max_segment_length=PARTIAL_SNAP_MAX_SEGMENT_LENGTH,
        threshold_exclusion_area=0,
        threshold_exclusion_percentage=0,
        threshold_inclusion_percentage=100,
        buffer_multiplication_factor=BUFFER_MULTIPLICATION_FACTOR,
        threshold_circle_ratio=0.98,
        correction_distance=0.01,
        diff_metric=DIFF_METRIC,
        mitre_limit=10,
        area_limit=None,
        max_workers=None,
    ):
        """
        Initializes the Aligner object

        Args:
            feedback (object, optional): Feedback object that can be added to show
                feedback in QGIS. Defaults to None.
            relevant_distance (int, optional): The relevant distance (in meters) for
                processing. Defaults to 1.
            relevant_distances ([],optional): Relevant distances (in meters) for
                processing
            od_strategy (int, optional): The strategy to determine how to handle
                information outside the reference polygons (Open Domain)
                (default: SNAP_FULL_AREA_ALL_SIDE)
            threshold_overlap_percentage (int, optional): Threshold (%) to determine
                from which overlapping-percentage a reference-polygon has to be included
                when there aren't relevant intersections or relevant differences
                (default 50%).
                When setting this parameter to '-1' the original border for will be returned for cases where nor relevant intersections and relevant differences are found
            crs (str, optional): Coordinate Reference System (CRS) of the data.
                (default EPSG:31370)
            multi_as_single_modus (boolean, optional): Modus to handle multipolygons (Default=True):
                True: input-multipolygons will be split-up into single polygons and handled by the algorithm. After executing the algorithm, the results are merged together.
                False: Multipolygons are directly processed by the algorithm
            threshold_exclusion_percentage (int, optional): Percentage for excluding candidate reference-polygons when overlap(%) is smaller than the threshold(Default=0)
            threshold_exclusion_area (int, optional):Area in m² for excluding candidate reference-polygons when overlap(m²) is smaller than the threshold (Default=0)
            buffer_multiplication_factor (float, optional): Multiplication factor, used to buffer the thematic objects searching for reference borders (buffer= buffer_multiplication_factor*relevant_distance)(Default=1.01)
            threshold_circle_ratio (float, optional): Threshold-value to exclude circles getting processed (perfect circle = 1) based on POLSBY-POPPER algorithm(Default=0.98)
            correction_distance (float, optional): Distance used in a pos_neg_buffer to remove slivers (technical correction) (Default= 0.01 = 1cm )
            mitre_limit (int, optional):buffer-parameter - The mitre ratio is the ratio of the distance from the corner to the end of the mitred offset corner.
                When two line segments meet at a sharp angle, a miter join will extend far beyond the original geometry. (and in the extreme case will be infinitely far.) To prevent unreasonable geometry, the mitre limit allows controlling the maximum length of the join corner.
                Corners with a ratio which exceed the limit will be beveled(Default=10)
            area_limit (int, optional): Maximum area for processing. (default 100000)
            max_workers (int, optional): Amount of workers that is used in ThreadPoolExecutor (for parallel execution) when processing objects for multiple relevant distances. (default None). If set to -1, no parallel exececution is used.

        """
        self.logger = Logger(feedback)
        if relevant_distances is None and relevant_distance is not None:
            relevant_distances = [relevant_distance]
        self.relevant_distances = relevant_distances
        self.od_strategy = od_strategy
        self.threshold_overlap_percentage = threshold_overlap_percentage
        # Area in m² for excluding candidate reference-polygons when overlap(m²) is smaller than the
        # threshold
        self.threshold_exclusion_area = threshold_exclusion_area
        # Percentage for excluding candidate reference-polygons when overlap(%) is smaller than the
        # threshold
        self.threshold_exclusion_percentage = threshold_exclusion_percentage
        self.threshold_inclusion_percentage = threshold_inclusion_percentage
        self.area_limit = area_limit
        self.max_workers = max_workers
        # Multiplication-factor used in OD-strategy 2 (SNAP-BOTH SIDED) when calculating
        # OD-area to take into account
        self.buffer_multiplication_factor = buffer_multiplication_factor
        # Threshold-value to exclude circles getting processed (perfect circle = 1) based on
        # POLSBY-POPPER algorithm
        self.threshold_circle_ratio = threshold_circle_ratio
        # Distance used in a pos_neg_buffer to remove slivers (technical correction)
        self.correction_distance = correction_distance
        # Buffer parameters:
        # Distance to limit a buffered corner (MITER-join-style parameter)
        #   Explanation and examples:
        #   https://shapely.readthedocs.io/en/stable/reference/shapely.buffer.html
        #   https://postgis.net/docs/ST_Buffer.html
        self.mitre_limit = mitre_limit
        # quad_segments = 8 (by default in shapely)

        # PROCESSING DEFAULTS
        # thematic
        # name of the identifier-field of the thematic data (id has to be unique)
        self.name_thematic_id = ID_THEME_FIELD_NAME
        # dictionary to store all thematic geometries to handle
        self.dict_thematic: dict[any, BaseGeometry] = {}
        # dictionary to store properties of the reference-features (optional)
        self.dict_thematic_properties: dict[any, dict] = {}
        # Dict to store source-information of the thematic dictionary
        self.dict_thematic_source: dict[any, str] = {}
        # dictionary to store all unioned thematic geometries
        self.thematic_union = None

        # reference

        # name of the identifier-field of the reference data (id has to be unique,f.e
        # CAPAKEY for GRB-parcels)
        self.name_reference_id = ID_REFERENCE_FIELD_NAME
        # dictionary to store all reference geometries
        self.dict_reference: dict[any, BaseGeometry] = {}
        # dictionary to store properties of the reference-features (optional)
        self.dict_reference_properties: dict[any, dict] = {}
        # Dict to store source-information of the reference dictionary
        self.dict_reference_source: dict[any, str] = {}
        # to save a unioned geometry of all reference polygons; needed for calculation
        # in most OD-strategies
        self.reference_union = None

        # to save a the reference_elements (points and lines) that form the reference borders; needed for networkx-calculation
        # in most OD-strategies
        self.reference_elements = None

        # results

        # output-dictionaries (all results of process()), grouped by theme_id and relevant_distance
        self.dict_processresults: dict[any, dict[float, ProcessResult]] = {}
        # dictionary with the 'predicted' results, grouped by theme_id and relevant_distance
        self.dict_predictions: dict[any, dict[float, ProcessResult]] = {}
        # dictionary with the 'evaluated predicted' results, grouped by theme_id and relevant_distance
        self.dict_evaluated_predictions: dict[any, dict[float, ProcessResult]] = {}
        # dictionary with the 'evaluated predicted' properties, grouped by theme_id and relevant_distance
        self.dict_evaluated_predictions_properties: dict[any, dict[float, {}]] = {}

        # Coordinate reference system
        # thematic geometries and reference geometries are assumed to be in the same CRS
        # before loading into the Aligner. No CRS-transformation will be performed.
        # When loading data, CRS is expected to be a projected CRS with units in 'meter
        # (m)'.
        # Default EPSG:31370 (Lambert72), alternative: EPSG:3812 (Lambert2008)
        self.CRS = crs
        # this parameter is used to treat multipolygon as single polygons. So polygons
        # with ID splitter are separately evaluated and merged on result.
        self.multi_as_single_modus = multi_as_single_modus
        self.preserve_topology = preserve_topology
        self.snap_strategy = snap_strategy
        self.snap_max_segment_length = snap_max_segment_length
        self.partial_snapping = partial_snapping
        self.partial_snap_strategy = partial_snap_strategy
        self.partial_snap_max_segment_length = partial_snap_max_segment_length
        self.diff_metric = diff_metric
        self.logger.feedback_info("Aligner initialized")

    ##########LOADERS##########################
    ###########################################

    def load_thematic_data(self, loader: Loader):
        """
        Loads the thematic features into the aligner
        :param loader:
        :return:
        """
        self.dict_thematic, self.dict_thematic_properties, self.dict_thematic_source = (
            loader.load_data()
        )

        self.thematic_union = None
        # TODO reset all aligner variables too? fe dict_processresults etc

    def load_reference_data(self, loader: Loader):
        """
        Loads the reference features into the aligner, and prepares the reference-data for processing
        :param loader:
        :return:
        """
        (
            self.dict_reference,
            self.dict_reference_properties,
            self.dict_reference_source,
        ) = loader.load_data()
        self._prepare_reference_data()
        # TODO reset all aligner variables too? fe dict_processresults etc

    ##########PROCESSORS#######################
    ###########################################

    def _process_geometry_by_snap(
        self,
        input_geometry: BaseGeometry,
        ref_intersections_geoms,
        relevant_distance,
        snap_strategy,
        snap_max_segment_length,
    ) -> ProcessResult:
        snapped = []

        # Open Domain
        virtual_reference = Polygon()
        if self.od_strategy != OpenDomainStrategy.EXCLUDE:
            virtual_reference = self._create_virtual_reference(
                input_geometry, relevant_distance, False
            )

        if self.od_strategy == OpenDomainStrategy.EXCLUDE:
            pass
        elif self.od_strategy == OpenDomainStrategy.AS_IS:
            geom_od = safe_intersection(input_geometry, virtual_reference)
            snapped.append(geom_od)
        else:
            ref_intersections_geoms.append(virtual_reference)

        ref_geometrycollection = GeometryCollection(ref_intersections_geoms)

        snapped_geom = snap_geometry_to_reference(
            input_geometry,
            ref_geometrycollection,
            snap_strategy,
            snap_max_segment_length,
            relevant_distance,
        )
        snapped.append(snapped_geom)
        # union all
        geom_preresult = safe_unary_union(snapped)
        result_dict = self._postprocess_preresult(
            geom_preresult,
            input_geometry,
            GeometryCollection(),
            GeometryCollection(),
            relevant_distance,
        )
        return result_dict

    def _process_geometry_by_network(
        self,
        input_geometry: BaseGeometry,
        relevant_distance: float = 1,
        snap_strategy=SnapStrategy.NO_PREFERENCE,
    ) -> ProcessResult:

        input_geometry = to_multi(input_geometry)
        # determine the reference_elements
        input_geometry_buffered = buffer_pos(
            input_geometry, relevant_distance * self.buffer_multiplication_factor
        )
        # Former method to get all the reference elements in the surrounding area, replaced by optimized method below
        # ref_intersections = self.reference_items.take(
        #     self.reference_tree.query(input_geometry_buffered)
        # ).tolist()
        # ref_intersections_geoms = [
        #     self.dict_reference[key] for key in ref_intersections
        # ]
        # reference = safe_unary_union(
        #     extract_points_lines_from_geometry(
        #         GeometryCollection(ref_intersections_geoms)
        #     )
        # )
        # Optimized method to get all the reference elements in the surrounding area (reference elements only calculated once for the full aligner
        reference = safe_unary_union(
            safe_intersection(self._get_reference_elements(), input_geometry_buffered)
        )

        geom_processed_list = []
        if isinstance(input_geometry, (MultiPolygon)):
            geom_processed_list = []
            for polygon in input_geometry.geoms:
                exterior = polygon.exterior
                interiors = polygon.interiors
                exterior_processed = self._process_by_network(
                    exterior,
                    reference,
                    relevant_distance,
                    snap_strategy=snap_strategy,
                    close_output=True,
                )
                interiors_processed = []
                for i in interiors:
                    i_processed = self._process_by_network(
                        i,
                        reference,
                        relevant_distance,
                        snap_strategy=snap_strategy,
                        close_output=True,
                    )
                    interiors_processed.append(i_processed)
                geom_processed = Polygon(exterior_processed, interiors_processed)
                geom_processed_list.append(geom_processed)

        else:
            for geom in input_geometry.geoms:
                geom_processed = self._process_by_network(
                    geom,
                    reference,
                    relevant_distance,
                    snap_strategy=snap_strategy,
                    close_output=False,
                )
                geom_processed_list.append(geom_processed)
        geom_processed = safe_unary_union(geom_processed_list)

        return self._postprocess_preresult(
            geom_processed,
            input_geometry,
            GeometryCollection(),
            GeometryCollection(),
            relevant_distance,
        )

    def _process_by_network(
        self,
        geom_to_process,
        reference,
        relevant_distance,
        snap_strategy=SnapStrategy.NO_PREFERENCE,
        close_output=False,
    ):
        """

        :param geom_to_process: only single geometries are allowed (No Multi)
        :param reference:
        :param relevant_distance:
        :param snap_strategy:
        :param close_output:
        :return:
        """
        buffer_distance = relevant_distance / 2
        geom_to_process_buffered = buffer_pos(geom_to_process, buffer_distance)
        reference_buffered = buffer_pos(reference, buffer_distance)

        overlap = safe_intersection(geom_to_process_buffered, reference_buffered)
        if overlap.is_empty:
            return geom_to_process  # processresult original
        overlap_buffered = buffer_pos(overlap, buffer_distance)

        reference_intersection = safe_intersection(reference, overlap_buffered)
        # reference_intersection=set_precision(reference_intersection,precision)
        reference_intersection = safe_unary_union(reference_intersection)
        reference_intersection = to_multi(reference_intersection)
        if reference_intersection.is_empty:
            return geom_to_process  # processresult original
        # thematic_intersection = safe_intersection(geom_to_process, overlap_buffered)
        # thematic_intersection = line_merge(safe_unary_union(thematic_intersection))
        # thematic_intersection = to_multi(thematic_intersection)

        thematic_difference = safe_difference(geom_to_process, overlap_buffered)
        # thematic_difference = set_precision(thematic_difference,precision)
        thematic_difference = safe_unary_union(thematic_difference)
        thematic_difference = to_multi(thematic_difference)

        # Calculate vertices of the reference_intersection
        if snap_strategy != SnapStrategy.NO_PREFERENCE:
            reference_coords = MultiPoint(list(get_coords_from_geometry(reference)))
            reference_coords_intersection = to_multi(
                safe_intersection(reference_coords, overlap_buffered)
            )
        else:
            reference_coords_intersection = None
        # Logic that manages the composition of the lines (connections etc) based on the geom-inputtype
        # POINT
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
        # LINESTRING OR POLYGON
        else:
            geom_processed = self._get_processed_network_path(
                geom_to_process,
                reference_intersection,
                reference_coords_intersection,
                thematic_difference,
                snap_strategy,
                relevant_distance,
            )

            if (
                close_output
                and geom_processed is not None
                and not geom_processed.is_ring
            ):
                closed_coords = list(geom_processed.coords) + [geom_processed.coords[0]]
                geom_processed = LineString(closed_coords)
        return make_valid(geom_processed)

    def _get_processed_network_path(
        self,
        geom_to_process,
        reference_intersection,
        reference_coords_intersection,
        thematic_difference,
        snap_strategy,
        relevant_distance,
    ):
        thematic_points = self._get_thematic_points(
            geom_to_process, reference_intersection
        )
        segments = []
        segments.extend(list(reference_intersection.geoms))
        segments.extend(list(thematic_difference.geoms))

        # add extra segments (connectionlines between theme and reference)
        extra_segments = []
        for geom in thematic_difference.geoms:
            p_start = Point(geom.coords[0])
            p_end = Point(geom.coords[-1])
            connection_line_start = self._get_connection_line(
                thematic_points,
                thematic_difference,
                p_start,
                reference_intersection,
                reference_coords_intersection,
                relevant_distance,
            )
            connection_line_end = self._get_connection_line(
                thematic_points,
                thematic_difference,
                p_end,
                reference_intersection,
                reference_coords_intersection,
                relevant_distance,
            )
            extra_segments.append(connection_line_start)
            extra_segments.append(connection_line_end)
        segments.extend(extra_segments)

        # add extra segments (connection lines between reference_intersections)
        # TODO, when passing OpenDomain, there is not always connection between the reference-parts. How to solve this?
        extra_segments_ref_intersections = shortest_connections_between_geometries(
            reference_intersection
        )
        # extra_segments_ref_intersections = get_connection_lines_to_nearest(
        #     reference_intersection
        # )
        segments.extend(
            extra_segments_ref_intersections
        )  # removed as we first going to filter these lines

        # Filter out lines that are not fully in relevant distance as these are no valid solution-paths
        # Mostly these lines will already be in the distance as both start en endpoint are in range (but not always fully)
        # geom_to_process_buffered = buffer_pos(geom_to_process, relevant_distance*1.01)
        # extra_geomcollection_ref_intersections = safe_intersection(
        #     GeometryCollection(extra_segments_ref_intersections),
        #     geom_to_process_buffered,
        # )
        # segments.append(extra_geomcollection_ref_intersections)

        network = prepare_network(segments)
        geom_processed = find_best_path_in_network(
            geom_to_process, network, snap_strategy, relevant_distance
        )
        # if geom_processed is None:
        #     # add original so a connected path will be found
        #     segments.append(geom_to_process)
        #     network = prepare_network(segments)
        #     geom_processed = find_longest_path_in_network(
        #         geom_to_process, network, snap_strategy, relevant_distance
        #     )

        return geom_processed

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
        geom_to_process_segmentized = segmentize(
            geom_to_process_line, self.partial_snap_max_segment_length
        )
        # Split the line at all intersection points with the MultiLineString
        splitter = safe_unary_union(reference_intersection)
        try:
            geom_to_process_splitted = split(geom_to_process_segmentized, splitter)
        except GeometryTypeError:
            geom_to_process_splitted = geom_to_process_segmentized
        thematic_points = MultiPoint(
            list(get_coords_from_geometry(geom_to_process_splitted))
        )
        return thematic_points

    def _get_connection_line(
        self,
        thematic_points,
        thematic_difference,
        point,
        reference_intersection,
        reference_coords_intersection,
        relevant_distance,
    ):
        # factor = 1.001

        closest_points = heapq.nsmallest(
            2, thematic_points.geoms, key=lambda p: point.distance(p)
        )
        points = [closest_points[0], point, closest_points[-1]]
        line_theme = LineString(points)
        dist_1 = closest_points[0].distance(thematic_difference)
        dist_2 = closest_points[-1].distance(thematic_difference)
        if dist_1 > dist_2:
            line_theme_furthest_point = closest_points[0]
        elif dist_1 < dist_2:
            line_theme_furthest_point = closest_points[-1]
        else:
            line_theme_furthest_point = point

        if (
            not reference_coords_intersection is None
            and not reference_coords_intersection.is_empty
        ):
            line_ref = shortest_line(
                line_theme_furthest_point, reference_coords_intersection
            )

            if line_ref.length > relevant_distance * 1.5:
                line_ref = shortest_line(
                    line_theme_furthest_point, reference_intersection
                )
        else:
            line_ref = shortest_line(line_theme_furthest_point, reference_intersection)

        # connection_line = scale(connection_line, factor, factor)
        # line_theme = scale(line_theme, factor, factor)
        # line_ref= scale(line_ref, factor, factor)
        connection_line = safe_unary_union([line_theme, line_ref])
        # To scale or not to scale, that's the question.
        # When vertices are used it is possibly not necessary to scale because we use the vertices of the segmentized input_geometry, so no problem with floating point-intersections.
        # When we do not use vertices it could be necessary (due to floating point error) to make sure lines are intersecting so they are split on these intersecting points

        if (
            round(connection_line.length, RELEVANT_DISTANCE_DECIMALS)
            > relevant_distance * self.partial_snap_max_segment_length * 2
            # * factor
            # * factor
            # * 4  # There could be a better way to exclude invalid connection-lines?
        ):
            return LineString()

        return connection_line

    def process_geometry(
        self,
        input_geometry: BaseGeometry,
        relevant_distance: float = 1,
        od_strategy=OpenDomainStrategy.SNAP_ALL_SIDE,
        threshold_overlap_percentage=50,
    ) -> ProcessResult:
        """
        method to align a geometry to the reference layer

        Args:
            input_geometry (BaseGeometry): The input geometric object.
            relevant_distance: The relevant distance (in meters) for processing
            od_strategy (int, optional): The strategy to determine how to handle
                information outside the reference polygons (Open Domain)
                (default: SNAP_FULL_AREA_ALL_SIDE)
            threshold_overlap_percentage (int, optional): Threshold (%) to determine
                from which overlapping-percentage a reference-polygon has to be included
                when there aren't relevant intersections or relevant differences
                (default 50%).
                When setting this parameter to '-1' the original border for will be returned for cases where nor relevant intersections and relevant differences are found


        Returns:
            ProcessResult : A dict containing the resulting geometries:

            *   result (BaseGeometry): The resulting output geometry
            *   result_diff (BaseGeometry): The resulting difference output geometry
            *   result_diff_plus (BaseGeometry): The resulting positive difference
                output geometry
            *   result_diff_min (BaseGeometry): The resulting negative difference output
                geometry
            *   relevant_intersection (BaseGeometry): The relevant_intersection
            *   relevant_difference (BaseGeometry): The relevant_difference
            *   remark (str): remarks collected when processing the geometry
        """
        # Processing based on thematic geom_type and reference_geom_type

        self.od_strategy = od_strategy
        self.threshold_overlap_percentage = threshold_overlap_percentage

        if isinstance(input_geometry, GeometryCollection):
            raise ValueError(
                "GeometryCollection as input is not supported. Please use the individual geometries from the GeometryCollection as input."
            )
        elif isinstance(input_geometry, (Polygon, MultiPolygon)):
            # Processing thematic polygons
            if self.area_limit and input_geometry.area > self.area_limit:
                message = f"The input polygon is too large to process: input area {str(input_geometry.area)} m², limit area: {str(self.area_limit)} m²."
                raise ValueError(message)
            self.logger.feedback_debug("process geometry")

            # CALCULATE INNER and OUTER INPUT GEOMETRY for performance optimisation on big geometries
            # combine all parts of the input geometry to one polygon
            input_geometry_inner, input_geometry_outer = _calculate_inner_outer(
                input_geometry, relevant_distance
            )
            # get a list of all ref_ids that are intersecting the thematic geometry; we take it bigger because we want to check if there are also reference geometries on a relevant distance.
            input_geometry_outer_buffered = buffer_pos(
                input_geometry_outer,
                relevant_distance * self.buffer_multiplication_factor,
            )
            ref_intersections = self.reference_items.take(
                self.reference_tree.query(input_geometry_outer_buffered)
            ).tolist()

            ref_intersections_geoms = []
            all_polygons = True
            for key_ref in ref_intersections:
                ref_geom = self.dict_reference[key_ref]
                ref_intersections_geoms.append(ref_geom)
                if not isinstance(ref_geom, (Polygon, MultiPolygon)):
                    all_polygons = False

            if all_polygons:
                return self._process_geometry_by_brdr(
                    input_geometry,
                    input_geometry_outer,
                    input_geometry_inner,
                    ref_intersections_geoms,
                    relevant_distance=relevant_distance,
                )
        return self._process_geometry_by_network(
            input_geometry,
            relevant_distance,
            snap_strategy=self.partial_snap_strategy,
        )

    def _process_geometry_by_brdr(
        self,
        input_geometry,
        input_geometry_outer,
        input_geometry_inner,
        ref_intersections_geoms,
        relevant_distance,
    ):
        buffer_distance = relevant_distance / 2

        # array with all relevant parts of a thematic geometry; initial empty Polygon
        (
            preresult,
            relevant_intersection_array,
            relevant_diff_array,
        ) = self._calculate_intersection_between_geometry_and_od(
            input_geometry_outer, input_geometry_inner, relevant_distance
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
            ) = _calculate_geom_by_intersection_and_reference(
                geom_intersection,
                geom_reference,
                input_geometry_inner,
                False,
                buffer_distance,
                self.threshold_overlap_percentage,
                self.threshold_exclusion_percentage,
                self.threshold_exclusion_area,
                self.threshold_inclusion_percentage,
                self.mitre_limit,
                self.partial_snapping,
                self.partial_snap_strategy,
                self.partial_snap_max_segment_length,
            )
            self.logger.feedback_debug("intersection calculated")

            preresult = self._add_multi_polygons_from_geom_to_array(geom, preresult)
            relevant_intersection_array = self._add_multi_polygons_from_geom_to_array(
                relevant_intersection, relevant_intersection_array
            )
            relevant_diff_array = self._add_multi_polygons_from_geom_to_array(
                relevant_diff, relevant_diff_array
            )

        # UNION INTERMEDIATE LAYERS
        relevant_intersection = safe_unary_union(relevant_intersection_array)
        if relevant_intersection is None or relevant_intersection.is_empty:
            relevant_intersection = Polygon()
        relevant_diff = safe_unary_union(relevant_diff_array)
        if relevant_diff is None or relevant_diff.is_empty:
            relevant_diff = Polygon()

        # Add inner input geometry to preresult
        preresult.append(input_geometry_inner)
        geom_preresult = safe_unary_union(preresult)

        # POSTPROCESSING
        result_dict = self._postprocess_preresult(
            geom_preresult,
            input_geometry,
            relevant_intersection,
            relevant_diff,
            relevant_distance,
        )

        return result_dict

    def process(
        self,
        dict_thematic=None,
        relevant_distances: Iterable[float] = None,
        relevant_distance=1,
        od_strategy=OpenDomainStrategy.SNAP_ALL_SIDE,
        threshold_overlap_percentage=50,
    ) -> dict[any, dict[float, ProcessResult]]:
        """
        Calculates the resulting dictionaries for thematic data based on a series of
            relevant distances.

        Args:
            relevant_distances (Iterable[float]): A series of relevant distances
                (in meters) to process
            od_strategy (int, optional): The strategy to determine how to handle
                information outside the reference polygons (Open Domain)
                (default: SNAP_FULL_AREA_ALL_SIDE)
            threshold_overlap_percentage (int, optional): Threshold (%) to determine
                from which overlapping-percentage a reference-polygon has to be included
                when there aren't relevant intersections or relevant differences
                (default 50%).
                When setting this parameter to '-1' the original border for will be returned for cases where nor relevant intersections and relevant differences are found


        Returns:
            dict: A dictionary, for every thematic ID a dictionary with the results for all distances

                {
                    'theme_id_1': {0: (ProcessResult), 0.1:
                        (ProcessResult), ...},
                    'theme_id_2': {0: (ProcessResult), 0.1:
                        (ProcessResult), ...},
                    ...
                }
        """
        if relevant_distances is None:
            relevant_distances = [relevant_distance]
        self.relevant_distances = relevant_distances
        self.od_strategy = od_strategy
        self.threshold_overlap_percentage = threshold_overlap_percentage
        self.logger.feedback_debug("Process series" + str(self.relevant_distances))
        dict_series = {}
        dict_series_queue = {}
        futures = []
        dict_thematic_to_process = dict_thematic
        if dict_thematic_to_process is None:
            dict_thematic_to_process = self.dict_thematic
        dict_multi_as_single = {}

        if self.multi_as_single_modus:
            dict_thematic_to_process, dict_multi_as_single = multi_to_singles(
                dict_thematic_to_process
            )

        if self.preserve_topology:
            # self.max_workers =-1
            # self.logger.feedback_info("max_workers set to -1 when using 'preserve_topology'")
            dict_thematic_to_process, topo_thematic = generate_topo(
                dict_thematic_to_process
            )

        if self.max_workers != -1:
            with ThreadPoolExecutor(
                max_workers=self.max_workers
            ) as executor:  # max_workers=5
                for key, geometry in dict_thematic_to_process.items():
                    self.logger.feedback_info(
                        f"thematic id {str(key)} processed with relevant distances (m) [{str(self.relevant_distances)}]"
                    )
                    dict_series[key] = {}
                    dict_series_queue[key] = {}
                    for relevant_distance in self.relevant_distances:
                        try:
                            future = executor.submit(
                                self.process_geometry,
                                geometry,
                                relevant_distance,
                                od_strategy,
                                threshold_overlap_percentage,
                            )
                            futures.append(future)
                            dict_series_queue[key][relevant_distance] = future
                        except ValueError as e:
                            self.logger.feedback_warning(
                                "error for"
                                + f"thematic id {str(key)} processed with relevant distances (m) [{str(self.relevant_distances)}]"
                            )
                            dict_series_queue[key][relevant_distance] = None
                            self.logger.feedback_warning(str(e))
            self.logger.feedback_debug("waiting all started RD calculations")
            wait(futures)
            for id_theme, dict_dist in dict_series_queue.items():
                for relevant_distance, future in dict_dist.items():
                    dict_series[id_theme][relevant_distance] = future.result()
        else:
            for key, geometry in dict_thematic_to_process.items():
                self.logger.feedback_info(
                    f"thematic id {str(key)} processed with relevant distances (m) [{str(self.relevant_distances)}]"
                )
                dict_series[key] = {}
                for relevant_distance in self.relevant_distances:
                    try:
                        processed_result = self.process_geometry(
                            geometry,
                            relevant_distance,
                            od_strategy,
                            threshold_overlap_percentage,
                        )
                    except ValueError as e:
                        self.logger.feedback_warning(str(e))
                        processed_result = None

                    dict_series[key][relevant_distance] = processed_result
        if self.preserve_topology:
            dict_series = dissolve_topo(
                dict_series,
                self.dict_thematic,
                dict_thematic_to_process,
                topo_thematic,
                relevant_distances,
            )

        if self.multi_as_single_modus:
            dict_series = merge_process_results(dict_series, dict_multi_as_single)

        # Check if geom changes from multi to single or vice versa
        for theme_id, dict_dist_results in dict_series.items():
            original_geometry = self.dict_thematic[theme_id]
            try:
                original_geometry_length = len(original_geometry.geoms)
            except:
                original_geometry_length = 1
            for relevant_distance, process_result in dict_dist_results.items():
                resulting_geom = process_result["result"]
                try:
                    resulting_geometry_length = len(resulting_geom.geoms)
                except:
                    resulting_geometry_length = 1
                if original_geometry_length != resulting_geometry_length:
                    msg = "Difference in amount of geometries"
                    self.logger.feedback_debug(msg)
                    process_result["remark"] = process_result["remark"] + " | " + msg

        self.logger.feedback_info(
            "End of processing series: " + str(self.relevant_distances)
        )
        self.dict_processresults = dict_series

        return self.dict_processresults

    def predictor(
        self,
        dict_thematic=None,
        relevant_distances=[
            round(k, RELEVANT_DISTANCE_DECIMALS)
            for k in np.arange(0, 310, 10, dtype=int) / 100
        ],
        od_strategy=OpenDomainStrategy.SNAP_ALL_SIDE,
        threshold_overlap_percentage=50,
        diff_metric=None,
    ):
        """
        Predicts the 'most interesting' relevant distances for changes in thematic
        elements based on a distance series.

        This function analyzes a set of thematic geometries (`self.dict_thematic`) to
        identify potentially interesting distances where changes occur. It performs
        the following steps:

        1. **Process Distance Series:**
            - Calculates a series of results for different distances specified by
              `relevant_distances`.
            - This calculation might involve functions like `self.process_series`
              (implementation details likely depend on your specific code).

        2. **Calculate Difference Metrics:**
            - Analyzes the results from the distance series to compute difference
              metrics between thematic elements at each distance (using
              `diffs_from_dict_series`).

        3. **Identify Breakpoints and Zero-Streaks:**
            - For each thematic geometry, it identifies potential "breakpoints" where
              the difference metric changes sign (from positive to negative or vice
              versa).
            - It also identifies "zero-streaks" which are consecutive distances with a
              difference metric close to zero (potentially indicating minimal change).

        4. **Predict Interesting Distances:**
            - The function considers distances corresponding to breakpoints and
              zero-streaks as potentially interesting for further analysis.
            - These distances are stored in a dictionary (`dict_predictions`) with the
              thematic element key as the outer key.
            - Additionally, the corresponding results from the distance series for
              those distances are included.

        Args:
            dict_thematic: the dictionary with the thematic geometries to 'predict'. Default is None, so all thematic geometries inside the aligner will be processed.
            relevant_distances (np.ndarray, optional): A series of relevant distances
                (in meters) to process. : A NumPy array of distances to
              be analyzed.
            od_strategy (int, optional): The strategy to determine how to handle
                information outside the reference polygons (Open Domain)
                (default: SNAP_FULL_AREA_ALL_SIDE)
            threshold_overlap_percentage (int, optional): Threshold (%) to determine
                from which overlapping-percentage a reference-polygon has to be included
                when there aren't relevant intersections or relevant differences
                (default 50%).
                When setting this parameter to '-1' the original border for will be returned for cases where nor relevant intersections and relevant differences are found

            relevant_distances (np.ndarray, optional): A NumPy array of distances to
              be analyzed. Defaults to np.arange(0.1, 5.05, 0.1).
            od_strategy (OpenDomainStrategy, optional): A strategy for handling
              open data in the processing (implementation specific). Defaults to
             OpenDomainStrategy.SNAP_ALL_SIDE.
            threshold_overlap_percentage (int, optional): A percentage threshold for
              considering full overlap in the processing (implementation specific).
             Defaults to 50.
            diff_metric (enum, optional): A enum thjat determines the method how differences are measured to determine the 'predictions'

        Returns:
            dict_series: A dictionary containing the resultset for all relevant distances for each thematic element.
            dict_predictions: A dictionary containing predicted interesting distances for each
            thematic element.
                - Keys: Thematic element identifiers from `self.dict_thematic`.
                - Values: Dictionaries with the following structure for each
                   thematic element:
                    - Keys: Distances identified as interesting (breakpoints or
                    zero-streaks).
                    - Values: dicts containing results (likely specific to
                    your implementation) from the distance series for the
                    corresponding distance.
            diffs_dict: a dictionary with the differences for each relevant distance
        """

        if dict_thematic is None:
            dict_thematic = self.dict_thematic
        dict_predictions = defaultdict(dict)
        if od_strategy is None:
            od_strategy = OpenDomainStrategy.SNAP_ALL_SIDE
        if threshold_overlap_percentage is None:
            threshold_overlap_percentage = 50
        relevant_distances = list(relevant_distances)
        relevant_distances.append(round(0, RELEVANT_DISTANCE_DECIMALS))
        relevant_distances = list(set(relevant_distances))
        relevant_distances = sorted(relevant_distances)
        dict_processresults = self.process(
            dict_thematic=dict_thematic,
            relevant_distances=relevant_distances,
            od_strategy=od_strategy,
            threshold_overlap_percentage=threshold_overlap_percentage,
        )
        if diff_metric is None:
            diff_metric = self.diff_metric

        diffs_dict = self.get_diff_metrics(
            dict_processresults, dict_thematic, diff_metric=diff_metric
        )

        for theme_id, diffs in diffs_dict.items():
            if len(diffs) != len(relevant_distances):
                self.logger.feedback_warning(
                    f"Number of computed diffs for thematic element {theme_id} does "
                    f"not match the number of relevant distances."
                )
                continue
            diff_values = list(diffs.values())
            breakpoints, zero_streaks = get_breakpoints_zerostreak(
                relevant_distances, diff_values
            )
            self.logger.feedback_debug(str(theme_id))
            if len(zero_streaks) == 0:
                self.logger.feedback_debug(
                    "No zero-streaks found for: " + str(theme_id)
                )
            for zs in zero_streaks:
                dict_predictions[theme_id][zs[0]] = dict_processresults[theme_id][zs[0]]
                dict_predictions[theme_id][zs[0]][PREDICTION_SCORE] = zs[3]

        # Check if the predicted reldists are unique (and remove duplicated predictions
        dict_predictions_unique = defaultdict(dict)
        for theme_id, dist_results in dict_predictions.items():
            dict_predictions_unique[theme_id] = {}
            predicted_geoms_for_theme_id = []
            for rel_dist, processresults in dist_results.items():
                predicted_geom = processresults["result"]
                if not _equal_geom_in_array(
                    predicted_geom,
                    predicted_geoms_for_theme_id,
                    self.correction_distance,
                    self.mitre_limit,
                ) or predicted_geom.geom_type in (
                    "Point",
                    "MultiPoint",
                    "LineString",
                    "MultiLineString",
                ):
                    dict_predictions_unique[theme_id][rel_dist] = processresults
                    predicted_geoms_for_theme_id.append(processresults["result"])
                else:
                    self.logger.feedback_info(
                        f"Duplicate prediction found for key {theme_id} at distance {rel_dist}: Prediction excluded"
                    )
            for dist in dict_predictions_unique[theme_id].keys():
                dict_predictions_unique[theme_id][dist][PREDICTION_COUNT] = len(
                    predicted_geoms_for_theme_id
                )

        self.dict_predictions = dict_predictions_unique

        return (
            dict_processresults,
            self.dict_predictions,
            diffs_dict,
        )

    def evaluate(
        self,
        ids_to_evaluate=None,
        base_formula_field=FORMULA_FIELD_NAME,
        relevant_distances=[
            round(k, RELEVANT_DISTANCE_DECIMALS)
            for k in np.arange(0, 310, 10, dtype=int) / 100
        ],
        full_strategy=FullStrategy.NO_FULL,
        max_predictions=-1,
        multi_to_best_prediction=True,
    ):
        """
        Compares and evaluate input-geometries (with formula). Attributes are added to evaluate and decide if new
        proposals can be used
        ids_to_evaluate: list with all IDs to evaluate. all other IDs will be unchanged. If None (default), all self.dict_thematic will be evaluated.
        base_formula_field: name of the field where the base_formula is found in the data
        max_predictions: integer that indicates how many predictions are maximally returned. (-1 indicates all predictions are returned)
        relevant_distances: relevant distances to evaluate
        full_strategy: enum, decided which predictions are kept or prefered based on full-ness of the prediction
        multi_to_best_prediction (default True): Only usable in combination with max_predictions=1. If True (and max_predictions=1), the prediction with highest score will be taken.If False, the original geometry is returned.
        """
        if ids_to_evaluate is None:
            ids_to_evaluate = list(self.dict_thematic.keys())
        dict_affected = {}
        dict_unaffected = {}
        for id_theme, geom in self.dict_thematic.items():
            if id_theme in ids_to_evaluate:
                dict_affected[id_theme] = geom
            else:
                dict_unaffected[id_theme] = geom
        # AFFECTED
        dict_series, dict_affected_predictions, diffs = self.predictor(
            dict_thematic=dict_affected,
            relevant_distances=relevant_distances,
            od_strategy=self.od_strategy,
            threshold_overlap_percentage=self.threshold_overlap_percentage,
            diff_metric=self.diff_metric,
        )
        dict_predictions_evaluated = {}
        prop_dictionary = {}

        for theme_id in dict_affected.keys():
            dict_predictions_evaluated[theme_id] = {}
            prop_dictionary[theme_id] = {}
            if theme_id not in dict_affected_predictions.keys():
                # No predictions available
                relevant_distance = round(0, RELEVANT_DISTANCE_DECIMALS)
                props = self._evaluate(
                    id_theme=theme_id,
                    geom_predicted=dict_affected[theme_id],
                    base_formula_field=base_formula_field,
                )
                props[EVALUATION_FIELD_NAME] = Evaluation.TO_CHECK_NO_PREDICTION
                props[PREDICTION_COUNT] = 0
                props[PREDICTION_SCORE] = -1
                dict_predictions_evaluated[theme_id][relevant_distance] = {
                    "result": dict_affected[theme_id],
                    "remark": "no prediction, original returned",
                }
                prop_dictionary[theme_id][relevant_distance] = props
                continue

            # When there are predictions available
            dict_predictions_results = dict_affected_predictions[theme_id]
            scores = []
            distances = []
            predictions = []
            # Moet dit niet naar de predictor verhuizen?
            prediction_properties = []
            equality_found = False

            for dist in sorted(dict_predictions_results.keys()):
                if equality_found:
                    continue
                # prop_dictionary[theme_id][dist] = {}
                props = self._evaluate(
                    id_theme=theme_id,
                    geom_predicted=dict_predictions_results[dist]["result"],
                    base_formula_field=base_formula_field,
                )
                prediction_count = dict_affected_predictions[theme_id][dist][
                    PREDICTION_COUNT
                ]
                prediction_score = dict_affected_predictions[theme_id][dist][
                    PREDICTION_SCORE
                ]
                props[PREDICTION_COUNT] = prediction_count
                props[PREDICTION_SCORE] = prediction_score
                full = props[FULL_ACTUAL_FIELD_NAME]
                if full_strategy == FullStrategy.ONLY_FULL and not full:
                    continue
                if (
                    props[EVALUATION_FIELD_NAME] == Evaluation.TO_CHECK_NO_PREDICTION
                    and props[PREDICTION_COUNT] == 1
                ):
                    props[EVALUATION_FIELD_NAME] = Evaluation.PREDICTION_UNIQUE
                elif (
                    props[EVALUATION_FIELD_NAME] == Evaluation.TO_CHECK_NO_PREDICTION
                    and props[PREDICTION_COUNT] > 1
                ):
                    props[EVALUATION_FIELD_NAME] = Evaluation.TO_CHECK_PREDICTION_MULTI
                elif props[EVALUATION_FIELD_NAME] != Evaluation.TO_CHECK_NO_PREDICTION:
                    dict_predictions_evaluated[theme_id][dist] = (
                        dict_affected_predictions[theme_id][dist]
                    )
                    prop_dictionary[theme_id][dist] = props
                    equality_found = True
                    continue
                if full:
                    if full_strategy != FullStrategy.NO_FULL:
                        props[EVALUATION_FIELD_NAME] = (
                            Evaluation.TO_CHECK_PREDICTION_FULL
                        )
                        prediction_score = prediction_score + 50
                        if prediction_score > 100:
                            prediction_score = 100
                        props[PREDICTION_SCORE] = prediction_score
                        props[PREDICTION_COUNT] = prediction_count
                    else:
                        props[EVALUATION_FIELD_NAME] = (
                            Evaluation.TO_CHECK_PREDICTION_MULTI_FULL
                        )
                scores.append(props[PREDICTION_SCORE])
                distances.append(dist)
                predictions.append(dict_affected_predictions[theme_id][dist])
                prediction_properties.append(props)

            # get max amount of best-scoring predictions
            best_ix = sorted(range(len(scores)), reverse=True, key=lambda i: scores[i])
            len_best_ix = len(best_ix)

            # if there is only one prediction left,  evaluation is set to PREDICTION_UNIQUE_FULL
            if len_best_ix == 1:
                if (
                    FULL_ACTUAL_FIELD_NAME in prediction_properties[0]
                    and prediction_properties[0][FULL_ACTUAL_FIELD_NAME]
                ):
                    prediction_properties[0][
                        EVALUATION_FIELD_NAME
                    ] = Evaluation.PREDICTION_UNIQUE_FULL
                else:
                    prediction_properties[0][
                        EVALUATION_FIELD_NAME
                    ] = Evaluation.PREDICTION_UNIQUE

            # if there are multiple predictions, but we want only one and we ask for the original
            if (
                len_best_ix > 1
                and max_predictions == 1
                and not multi_to_best_prediction
            ):
                relevant_distance = round(0, RELEVANT_DISTANCE_DECIMALS)
                props[EVALUATION_FIELD_NAME] = Evaluation.TO_CHECK_ORIGINAL
                props[PREDICTION_SCORE] = -1
                dict_predictions_evaluated[theme_id][relevant_distance] = {
                    "result": dict_affected[theme_id],
                    "remark": "multiple predictions, original returned",
                }
                prop_dictionary[theme_id][relevant_distance] = props
                continue

            if max_predictions > 0 and len_best_ix > max_predictions:
                best_ix = best_ix[:max_predictions]
            if len(best_ix) > 0:
                for ix in best_ix:
                    distance = distances[ix]
                    prediction = predictions[ix]
                    props = prediction_properties[ix]
                    dict_predictions_evaluated[theme_id][distance] = prediction
                    prop_dictionary[theme_id][distance] = props
            else:
                # #when no evaluated predictions, the original is returned
                relevant_distance = round(0, RELEVANT_DISTANCE_DECIMALS)
                props = self._evaluate(
                    id_theme=theme_id,
                    geom_predicted=dict_affected[theme_id],
                    base_formula_field=base_formula_field,
                )
                props[EVALUATION_FIELD_NAME] = Evaluation.TO_CHECK_NO_PREDICTION
                props[PREDICTION_SCORE] = -1
                props[PREDICTION_COUNT] = 0
                dict_predictions_evaluated[theme_id][relevant_distance] = {
                    "result": dict_affected[theme_id],
                    "remark": "no prediction, original returned",
                }
                prop_dictionary[theme_id][relevant_distance] = props

        # UNAFFECTED
        relevant_distance = round(0, RELEVANT_DISTANCE_DECIMALS)
        for theme_id, geom in dict_unaffected.items():
            dict_predictions_evaluated[theme_id] = {}
            prop_dictionary[theme_id] = {relevant_distance: {}}
            props = self._evaluate(
                id_theme=theme_id,
                geom_predicted=geom,
                base_formula_field=base_formula_field,
            )
            props[EVALUATION_FIELD_NAME] = Evaluation.NO_CHANGE
            props[PREDICTION_SCORE] = -1
            dict_predictions_evaluated[theme_id][relevant_distance] = {"result": geom}
            prop_dictionary[theme_id][relevant_distance] = props
        self.dict_evaluated_predictions = dict_predictions_evaluated
        self.dict_evaluated_predictions_properties = prop_dictionary
        return dict_predictions_evaluated, prop_dictionary

    def get_brdr_formula(self, geometry: BaseGeometry, with_geom=False):
        """
        Calculates formula-related information based on the input geometry.

        Args:
            geometry (shapely.geometry object): The input geometry.
            with_geom (bool, optional): Whether to include geometry information in the
                output. Defaults to False.

        Returns:
            dict: A dictionary containing formula-related data:

            - "alignment_date": datetime.now().strftime(DATE_FORMAT),
            - "brdr_version": str(__version__),
            - "reference_source": self.dict_reference_source,
            - "full": True if the geometry exists out of all full reference-polygons, else False.
            - "area": Area of the geometry.
            - "reference_features": {
                array of all the reference features the geometry is composed of:
                    -   'full': True if the intersection is the same as the reference
                        geometry, else False.
                    -   'area': Area of the intersection or reference geometry.
                    -   'percentage': Percentage of intersection area relative to the
                        reference geometry.
                    -   'geometry': GeoJSON representation of the intersection (if
                        with_geom is True).},
            - "reference_od": Discription of the OD-part of the geometry (= not covered by reference-features),
        """
        dict_formula = {
            "alignment_date": datetime.now().strftime(DATE_FORMAT),
            "brdr_version": str(__version__),
            "reference_source": self.dict_reference_source,
            "full": True,
            "area": round(geometry.area, 2),
            "reference_features": {},
            "reference_od": None,
        }

        full_total = True
        last_version_date = None

        ref_intersections = self.reference_items.take(
            self.reference_tree.query(geometry)
        ).tolist()
        intersected = []
        for key_ref in ref_intersections:
            geom = None
            version_date = None
            geom_reference = self.dict_reference[key_ref]
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
                key_ref in self.dict_reference_properties
                and VERSION_DATE in self.dict_reference_properties[key_ref]
            ):
                str_version_date = self.dict_reference_properties[key_ref][VERSION_DATE]
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
        geom_od = buffer_pos(
            buffer_neg(
                safe_difference(geometry, safe_unary_union(intersected)),
                self.correction_distance,
                mitre_limit=self.mitre_limit,
            ),
            self.correction_distance,
            mitre_limit=self.mitre_limit,
        )
        if geom_od is not None:
            area_od = round(geom_od.area, 2)
            if area_od > 0:
                dict_formula["reference_od"] = {"area": area_od}
                if with_geom:
                    dict_formula["reference_od"]["geometry"] = to_geojson(geom_od)
        self.logger.feedback_debug(str(dict_formula))
        return dict_formula

    def get_diff_metrics(
        self,
        dict_processresults=None,
        dict_thematic=None,
        diff_metric=DiffMetric.CHANGES_AREA,
    ):
        """
        Calculates a dictionary containing difference metrics for thematic elements based on a distance series.

        Parameters:
        dict_series (dict): A dictionary where keys are thematic IDs and values are dictionaries mapping relative distances to ProcessResult objects.
        dict_thematic (dict): A dictionary where keys are thematic IDs and values are BaseGeometry objects representing the original geometries.
        diff_metric (DiffMetric, optional): The metric to use for calculating differences. Default is DiffMetric.CHANGES_AREA.

        Returns:
        dict: A dictionary where keys are thematic IDs and values are dictionaries mapping relative distances to calculated difference metrics.
        """
        if dict_processresults is None:
            dict_processresults = self.dict_processresults
        if dict_thematic is None:
            dict_thematic = self.dict_thematic
        return diffs_from_dict_processresults(
            dict_processresults=dict_processresults,
            dict_thematic=dict_thematic,
            reference_union=self._get_reference_union(),
            diff_metric=diff_metric,
        )

    ##########EXPORTERS########################
    ###########################################
    def get_results_as_geojson(
        self,
        resulttype=AlignerResultType.PROCESSRESULTS,
        formula=False,
        attributes=False,
    ):
        """
        get a geojson of a dictionary containing the resulting geometries for all
            'serial' relevant distances. The resulttype can be chosen.
        formula (boolean, Optional): The descriptive formula is added as an attribute to the result
        attributes (boolean, Optional): The original attributes/properties are added to the result
        """
        prop_dictionary = None
        if resulttype == AlignerResultType.PROCESSRESULTS:
            dict_series = self.dict_processresults
        elif resulttype == AlignerResultType.PREDICTIONS:
            dict_series = self.dict_predictions
        elif resulttype == AlignerResultType.EVALUATED_PREDICTIONS:
            dict_series = self.dict_evaluated_predictions
            prop_dictionary = self.dict_evaluated_predictions_properties
        else:
            raise (ValueError, "AlignerResultType unknown")
        if dict_series is None or dict_series == {}:
            self.logger.feedback_warning(
                "Empty results: No calculated results to export."
            )
            return {}
        if prop_dictionary is None:
            prop_dictionary = defaultdict(dict)

        for theme_id, results_dict in dict_series.items():
            for relevant_distance, process_results in results_dict.items():
                if relevant_distance not in prop_dictionary[theme_id]:
                    prop_dictionary[theme_id][relevant_distance] = {}
                if attributes and theme_id in self.dict_thematic_properties.keys():
                    for attr, value in self.dict_thematic_properties[theme_id].items():
                        prop_dictionary[theme_id][relevant_distance][attr] = value
                if (
                    formula
                ):  # and not (theme_id in prop_dictionary and relevant_distance in prop_dictionary[theme_id] and NEW_FORMULA_FIELD_NAME in prop_dictionary[theme_id][relevant_distance]):
                    result = process_results["result"]
                    formula = self.get_brdr_formula(result)
                    prop_dictionary[theme_id][relevant_distance][FORMULA_FIELD_NAME] = (
                        json.dumps(formula)
                    )
        return get_dict_geojsons_from_series_dict(
            dict_series,
            crs=self.CRS,
            id_field=self.name_thematic_id,
            series_prop_dict=prop_dictionary,
        )

    def get_input_as_geojson(self, inputtype=AlignerInputType.REFERENCE):
        """
        get a geojson of the input polygons (thematic or reference-polygons)
        """

        if inputtype == AlignerInputType.THEMATIC:
            dict_to_geojson = self.dict_thematic
            dict_properties = self.dict_thematic_properties
            property_id = self.name_thematic_id
        elif inputtype == AlignerInputType.REFERENCE:
            dict_to_geojson = self.dict_reference
            dict_properties = self.dict_reference_properties
            property_id = self.name_reference_id
        else:
            raise (ValueError, "AlignerInputType unknown")
        if dict_to_geojson is None or dict_to_geojson == {}:
            self.logger.feedback_warning("Empty input: No input to export.")
            return {}
        return geojson_from_dict(
            dict_to_geojson,
            self.CRS,
            property_id,
            prop_dict=dict_properties,
            geom_attributes=False,
        )

    def save_results(
        self, path, resulttype=AlignerResultType.PROCESSRESULTS, formula=True
    ):
        """
        Exports analysis results (as geojson) to path.

        This function exports 6 GeoJSON files containing the analysis results to the
        specified `path`.

        Args:
        path (str): The path to the directory where the GeoJSON files will be saved.
        formula (bool, optional): Whether to include formula-related information
            in the output. Defaults to True.

        Details of exported files:
        - result.geojson: Contains the original thematic data from `
          self.dict_result`.
        - result_diff.geojson: Contains the difference between the original
          and predicted data from `self.dict_result_diff`.
        - result_diff_plus.geojson: Contains results for areas that are
          added (increased area).
        - result_diff_min.geojson: Contains results for areas that are
          removed (decreased area).
        - result_relevant_intersection.geojson: Contains the areas with
          relevant intersection that has to be included in the result.
        - result_relevant_difference.geojson: Contains the areas with
          relevant difference that has to be excluded from the result.
        """

        fcs = self.get_results_as_geojson(
            formula=formula,
            resulttype=resulttype,
        )
        for name, fc in fcs.items():
            write_geojson(
                os.path.join(path, resulttype.value + "_" + name + ".geojson"), fc
            )

    def get_thematic_union(self):
        """
        returns a unary_unioned geometry from all the thematic geometries
        :return:
        """
        if self.thematic_union is None:
            self.thematic_union = safe_unary_union(list(self.dict_thematic.values()))
        return self.thematic_union

    def _prepare_reference_data(self):
        """
        Prepares reference data for spatial queries and analysis.
        It performs the following tasks:
        1. **Optimizes spatial queries:**
            - Creates a Spatial Relationship Tree (STRtree) using `STRtree` for
              efficient spatial queries against the reference data in
              `self.dict_reference`.
            - Converts the dictionary keys (reference identifiers) to a NumPy array
              for potential performance benefits in future operations.

        2. **Clears reference union:**
            - Sets `self.reference_union` to `None`. This variable stores the combined
              geometry of all reference data, and it's cleared here to indicate that
              it needs to be recalculated if requested later.

        Returns:
            None
        """
        # create an SRTree for performance optimisation
        self.logger.feedback_info(
            "length of reference_dict: " + str(len(self.dict_reference))
        )
        self.reference_tree = STRtree(list(self.dict_reference.values()))
        self.reference_items = np.array(list(self.dict_reference.keys()), dtype=object)
        # clear the reference_union, so it will be recalculated on request when needed
        self.reference_union = None
        return

    def _calculate_intersection_between_geometry_and_od(
        self, input_geometry, input_geometry_inner, relevant_distance
    ):
        """
        Calculates the intersecting parts between a thematic geometry and the Open Domain( domain, not coverd by reference-polygons)
        :param input_geometry:
        :param relevant_distance:
        :return:
        """
        # Calculate the intersection between thematic and Open Domain
        # buffer_distance = relevant_distance / 2
        relevant_intersection_array = []
        relevant_difference_array = []
        geom_thematic_od = Polygon()

        if self.od_strategy == OpenDomainStrategy.AS_IS:
            # All parts that are not covered by the reference layer are added to the
            #         resulting geometry AS IS
            self.logger.feedback_debug("OD-strategy AS IS")
            # all OD-parts wil be added AS IS
            geom_thematic_od = safe_difference(
                input_geometry, self._get_reference_union()
            )

        elif (
            self.od_strategy == OpenDomainStrategy.SNAP_INNER_SIDE
            or self.od_strategy == OpenDomainStrategy.EXCLUDE
        ):
            # integrates the entire inner area of the input geometry,
            # so Open Domain of the inner area is included in the result
            self.logger.feedback_debug("OD-strategy OD_SNAP_INNER_SIDE or EXCLUDE")
            geom_thematic_od = self._od_full_area(input_geometry, relevant_distance)

        elif self.od_strategy == OpenDomainStrategy.SNAP_ALL_SIDE:
            #  Inner & Outer-reference-boundaries are used.
            # integrates the entire inner area of the input geometry,
            self.logger.feedback_debug("OD-strategy OD-SNAP_ALL_SIDE")
            (
                geom_thematic_od,
                relevant_difference_array,
                relevant_intersection_array,
            ) = self._od_snap_all_side(
                input_geometry, input_geometry_inner, relevant_distance, outer=False
            )
            # This part calculates the full area
            geom_theme_od_min_clipped_plus_buffered_clipped = self._od_full_area(
                input_geometry, relevant_distance
            )
            # UNION of both elements
            geom_thematic_od = safe_union(
                geom_theme_od_min_clipped_plus_buffered_clipped, geom_thematic_od
            )
        elif self.od_strategy == OpenDomainStrategy.SNAP_PREFER_VERTICES:
            self.logger.feedback_debug("OD-strategy SNAP_PREFER_VERTICES")
            geom_thematic_od = self._od_snap(
                geometry=input_geometry,
                relevant_distance=relevant_distance,
                snap_strategy=SnapStrategy.PREFER_VERTICES,
            )

        elif self.od_strategy == OpenDomainStrategy.SNAP_NO_PREFERENCE:
            self.logger.feedback_debug("OD-strategy SNAP_NO_PREFERENCE")
            geom_thematic_od = self._od_snap(
                geometry=input_geometry,
                relevant_distance=relevant_distance,
                snap_strategy=SnapStrategy.NO_PREFERENCE,
            )

        elif self.od_strategy == OpenDomainStrategy.SNAP_ONLY_VERTICES:
            self.logger.feedback_debug("OD-strategy SNAP_ONLY_VERTICES")
            geom_thematic_od = self._od_snap(
                geometry=input_geometry,
                relevant_distance=relevant_distance,
                snap_strategy=SnapStrategy.ONLY_VERTICES,
            )

        # ADD THEMATIC_OD
        preresult = self._add_multi_polygons_from_geom_to_array(geom_thematic_od, [])
        return (
            preresult,
            relevant_intersection_array,
            relevant_difference_array,
        )

    def _od_snap(self, geometry, relevant_distance, snap_strategy):
        # all OD-parts wil be added AS IS
        geom_thematic_od = safe_difference(geometry, self._get_reference_union())
        if geom_thematic_od is None or geom_thematic_od.is_empty:
            return geom_thematic_od
        geom_thematic_od = to_multi(geom_thematic_od, geomtype="Polygon")
        out = []
        for p in geom_thematic_od.geoms:
            reference = safe_intersection(
                self._get_reference_union(),
                make_valid(
                    buffer_pos(
                        p,
                        self.buffer_multiplication_factor * relevant_distance,
                        mitre_limit=self.mitre_limit,
                    )
                ),
            )
            p_snapped = snap_geometry_to_reference(
                p,
                reference,
                max_segment_length=self.partial_snap_max_segment_length,
                snap_strategy=snap_strategy,
                tolerance=relevant_distance,
            )
            out.append(p_snapped)
        return safe_unary_union(out)

    def _od_full_area(self, geometry, relevant_distance):
        buffer_distance = relevant_distance / 2
        geom_theme_od = safe_difference(geometry, self._get_reference_union())
        geom_theme_min_buffered = buffer_neg(
            buffer_pos(
                buffer_neg(
                    geometry,
                    relevant_distance,
                    mitre_limit=self.mitre_limit,
                ),
                buffer_distance,
                mitre_limit=self.mitre_limit,
            ),
            buffer_distance,
            mitre_limit=self.mitre_limit,
        )
        geom_theme_od_clipped_min_buffered = safe_intersection(
            geom_theme_min_buffered, geom_theme_od
        )
        geom_theme_od_min_clipped_plus_buffered = buffer_pos(
            geom_theme_od_clipped_min_buffered,
            relevant_distance,
            mitre_limit=self.mitre_limit,
        )
        geom_theme_od_min_clipped_plus_buffered_clipped = safe_intersection(
            geom_theme_od_min_clipped_plus_buffered, geom_theme_od
        )
        geom_thematic_od = geom_theme_od_min_clipped_plus_buffered_clipped
        return geom_thematic_od

    def _od_snap_all_side(
        self, geometry, input_geometry_inner, relevant_distance, outer=False
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
            geometry, relevant_distance, outer
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
            ) = _calculate_geom_by_intersection_and_reference(
                geom_intersection,
                geom_reference,
                input_geometry_inner,
                True,
                buffer_distance,
                self.threshold_overlap_percentage,
                self.threshold_exclusion_percentage,
                self.threshold_exclusion_area,
                self.threshold_inclusion_percentage,
                self.mitre_limit,
                self.partial_snapping,
                self.partial_snap_strategy,
                self.partial_snap_max_segment_length,
            )

            relevant_intersection_array = self._add_multi_polygons_from_geom_to_array(
                geom_relevant_intersection, []
            )
            relevant_difference_array = self._add_multi_polygons_from_geom_to_array(
                geom_relevant_diff, []
            )
            # geom_thematic_od = safe_unary_union(geom_thematic_od_array)
        return geom_thematic_od, relevant_difference_array, relevant_intersection_array

    def _create_virtual_reference(self, geometry, relevant_distance, outer=False):
        """
        Functions that creates a 'virtual reference polygon' for areas that are not covered by reference-polygons.
        :param geometry:
        :param relevant_distance:
        :param outer: when outer is True, the outer boundary is used, inner is not used
        :return:
        """
        # TODO: remove outer, as it is not used?
        buffer_distance = relevant_distance / 2
        geom_thematic_buffered = make_valid(
            buffer_pos(
                geometry,
                (self.buffer_multiplication_factor * relevant_distance)
                + self.correction_distance,
                mitre_limit=self.mitre_limit,
            )
        )
        clip_ref_thematic_buffered = safe_intersection(
            self._get_reference_union(), geom_thematic_buffered
        )
        virtual_reference = safe_difference(
            geom_thematic_buffered, clip_ref_thematic_buffered
        )  # Both OD-parts are SNAPPED added
        if outer:  # when outer is True, the outer boundary is used, inner is not used
            geom_1 = safe_difference(geometry, virtual_reference)
            geom_2 = buffer_neg_pos(geom_1, buffer_distance)
            geom_3 = safe_intersection(geom_2, geometry)
            virtual_reference = safe_unary_union([geom_3, virtual_reference])
        return virtual_reference

    def _get_reference_union(self):
        """
        returns a unary_unioned geometry from all the reference geometries
        :return:
        """
        if self.reference_union is None:
            self.reference_union = safe_unary_union(list(self.dict_reference.values()))
        return self.reference_union

    def _get_reference_elements(self):
        """
        returns the points and lines from the reference geometries
        :return:
        """
        if self.reference_elements is None:
            self.reference_elements = extract_points_lines_from_geometry(
                GeometryCollection(list(self.dict_reference.values()))
            )
        return self.reference_elements

    def _postprocess_preresult(
        self,
        geom_preresult,
        geom_thematic,
        relevant_intersection,
        relevant_diff,
        relevant_distance,
    ) -> ProcessResult:
        """
        Postprocess the preresulting geometry with the following actions to create the final result
        *Corrections for areas that differ more than the relevant distance
        *slivers
        *Inner holes (donuts) /multipolygons
        *validity
        *Circles (Polsby-Popper)
        *Null/Empty-values

        Args:
            geom_preresult (BaseGeometry): The preresulting geometry to postprocess
            geom_thematic (BaseGeometry): The input geometry

        Returns:
            ProcessResult: A dictionary containing the resulting output geometries:

            *   result (BaseGeometry): The resulting output geometry
            *   result_diff (BaseGeometry): The resulting difference output geometry
            *   result_diff_plus (BaseGeometry): The resulting positive difference
                output geometry
            *   result_diff_min (BaseGeometry): The resulting negative difference output
                geometry
            *   remark (str): Remark when processing the geometry
        """
        remark = ""
        geom_thematic = make_valid(geom_thematic)
        if geom_preresult is None or geom_preresult.is_empty:
            # geom_preresult = geom_thematic
            geom_preresult = GeometryCollection()
            remark = "Empty geometry calculated: -->resulting geometry = empty geometry"
        if to_multi(geom_preresult).geom_type != to_multi(geom_thematic).geom_type:
            # geom_preresult = geom_thematic
            geom_preresult = GeometryCollection()
            remark = "Calculated geometry of different geomtype: -->resulting geometry = empty geometry"
        if geom_preresult.geom_type in [
            "Point",
            "MultiPoint",
            "LineString",
            "MultiLineString",
            "GeometryCollection",
        ]:
            result_diff_plus = safe_difference(
                geom_preresult, buffer_pos(geom_thematic, self.correction_distance)
            )
            result_diff_min = safe_difference(
                geom_thematic, buffer_pos(geom_preresult, self.correction_distance)
            )
            result_diff = safe_unary_union([result_diff_plus, result_diff_min])
            return _unary_union_result_dict(
                {
                    "result": geom_preresult,
                    "result_diff": result_diff,
                    "result_diff_plus": result_diff_plus,
                    "result_diff_min": result_diff_min,
                    "result_relevant_intersection": relevant_intersection,
                    "result_relevant_diff": relevant_diff,
                    "remark": remark,
                }
            )
        # Process array

        buffer_distance = relevant_distance / 2
        result = []
        geom_thematic_for_add_delete = geom_thematic

        if self.od_strategy == OpenDomainStrategy.EXCLUDE:
            geom_thematic_for_add_delete = safe_intersection(
                geom_thematic_for_add_delete, self._get_reference_union()
            )
            geom_preresult = safe_intersection(
                geom_preresult, self._get_reference_union()
            )

        if not (geom_thematic is None or geom_thematic.is_empty):
            # Correction for circles
            # calculate ratio to see if it is a circle, and keep the original geometry
            #  if a circle: (Polsby-Popper score)
            if (
                get_shape_index(geom_thematic.area, geom_thematic.length)
                > self.threshold_circle_ratio
            ):
                remark = "Circle detected: -->resulting geometry = original geometry"
                self.logger.feedback_debug(remark)
                return _unary_union_result_dict(
                    {"result": geom_thematic, "remark": remark}
                )

            # Correction for unchanged geometries
            # if safe_symmetric_difference(geom_preresult, geom_thematic).is_empty:
            if geometric_equality(
                geom_preresult,
                geom_thematic,
                correction_distance=self.correction_distance,
                mitre_limit=self.mitre_limit,
            ):
                remark = "Unchanged geometry: -->resulting geometry = original geometry"
                self.logger.feedback_debug(remark)
                return _unary_union_result_dict(
                    {"result": geom_thematic, "remark": remark}
                )

        # Corrections for areas that differ more than the relevant distance
        geom_thematic_dissolved = buffer_pos(
            buffer_neg(
                buffer_pos(
                    geom_preresult,
                    self.correction_distance,
                    mitre_limit=self.mitre_limit,
                ),
                2 * self.correction_distance,
                mitre_limit=self.mitre_limit,
            ),
            self.correction_distance,
            mitre_limit=self.mitre_limit,
        )
        # geom_symdiff = self._safe_symmetric_difference(geom_thematic,
        # geom_thematic_dissolved)
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
                    geom_diff_delete,
                    buffer_distance,
                    mitre_limit=self.mitre_limit,
                ),
            ),
        )
        geom_diff_removed_added = safe_union(
            geom_diff_removed,
            safe_intersection(
                geom_diff_add,
                buffer_neg_pos(
                    geom_diff_add,
                    buffer_distance,
                    mitre_limit=self.mitre_limit,
                ),
            ),
        )
        geom_thematic_preresult = buffer_pos(
            buffer_neg(
                buffer_pos(
                    geom_diff_removed_added,
                    self.correction_distance,
                    mitre_limit=self.mitre_limit,
                ),
                2 * self.correction_distance,
                mitre_limit=self.mitre_limit,
            ),
            self.correction_distance,
            mitre_limit=self.mitre_limit,
        )
        # Correction for Inner holes(donuts) / multipolygons
        # fill and remove gaps
        # TODO improvement: when relevant_distance very big (fe 100m) it could happen that parts of multipolygon-results will be removed unintentionally because part smaller than negative buffer
        geom_thematic_cleaned_holes = fill_and_remove_gaps(
            geom_thematic_preresult, buffer_distance
        )
        geom_thematic_result = buffer_pos(
            buffer_neg(
                buffer_pos(
                    geom_thematic_cleaned_holes,
                    self.correction_distance,
                    mitre_limit=self.mitre_limit,
                ),
                2 * self.correction_distance,
                mitre_limit=self.mitre_limit,
            ),
            self.correction_distance,
            mitre_limit=self.mitre_limit,
        )
        geom_thematic_result = make_valid(remove_repeated_points(geom_thematic_result))

        # Correction for empty preresults
        if geom_thematic_result.is_empty or geom_thematic_result is None:
            remark = "Calculated empty result: -->original geometry returned"
            self.logger.feedback_warning(remark)

            geom_thematic_result = geom_thematic
            # geom_thematic_result = Polygon() #If we return an empty geometry, the feature disappears, so we return the original geometry

        # group all initial multipolygons into a new resulting dictionary
        result.append(geom_thematic_result)

        # create all resulting geometries
        geom_thematic_result = safe_unary_union(result)

        # negative and positive buffer is added to the difference-calculations, to
        # remove 'very small' differences (smaller than the correction distance)
        geom_result_diff = buffer_pos(
            buffer_neg(
                safe_symmetric_difference(
                    geom_thematic_result,
                    geom_thematic,
                ),
                self.correction_distance,
                mitre_limit=self.mitre_limit,
            ),
            self.correction_distance,
            mitre_limit=self.mitre_limit,
        )
        geom_result_diff_plus = buffer_pos(
            buffer_neg(
                safe_difference(
                    geom_thematic_result,
                    geom_thematic,
                ),
                self.correction_distance,
                mitre_limit=self.mitre_limit,
            ),
            self.correction_distance,
            mitre_limit=self.mitre_limit,
        )
        geom_result_diff_min = buffer_pos(
            buffer_neg(
                safe_difference(
                    geom_thematic,
                    geom_thematic_result,
                ),
                self.correction_distance,
                mitre_limit=self.mitre_limit,
            ),
            self.correction_distance,
            mitre_limit=self.mitre_limit,
        )

        return _unary_union_result_dict(
            {
                "result": geom_thematic_result,
                "result_diff": geom_result_diff,
                "result_diff_plus": geom_result_diff_plus,
                "result_diff_min": geom_result_diff_min,
                "result_relevant_intersection": relevant_intersection,
                "result_relevant_diff": relevant_diff,
                "remark": remark,
            }
        )

    def _evaluate(
        self, id_theme, geom_predicted, base_formula_field=FORMULA_FIELD_NAME
    ):
        """
        function that evaluates a predicted geometry and returns a properties-dictionary
        """
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
        actual_formula = self.get_brdr_formula(geom_predicted)
        if actual_formula is None:
            properties[EVALUATION_FIELD_NAME] = Evaluation.TO_CHECK_NO_PREDICTION
            properties[FULL_ACTUAL_FIELD_NAME] = False
            return properties
        properties[FULL_ACTUAL_FIELD_NAME] = actual_formula["full"]
        properties[FORMULA_FIELD_NAME] = json.dumps(actual_formula)

        try:
            base_formula = json.loads(
                self.dict_thematic_properties[id_theme][base_formula_field]
            )
        except:
            base_formula = None

        if not is_brdr_formula(base_formula):
            properties[EVALUATION_FIELD_NAME] = Evaluation.TO_CHECK_NO_PREDICTION
            return properties
        properties[FULL_BASE_FIELD_NAME] = base_formula["full"]
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

        equal_reference_features = True
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
        # EVALUATION
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
        ):  # formula is the same,  both geometries are not 'full'
            properties[EVALUATION_FIELD_NAME] = Evaluation.EQUALITY_EQUAL_FORMULA_2
        elif (
            base_formula["full"] and actual_formula["full"] and od_alike
        ):  # formula not the same but geometries are full
            properties[EVALUATION_FIELD_NAME] = Evaluation.EQUALITY_FULL_3
        # elif base_formula["full"] == actual_formula["full"] and od_alike:
        #    properties[EVALUATION_FIELD_NAME] = Evaluation.EQUALITY_NON_FULL
        # TODO evaluate when not-full-parcels: compare all parcels?
        # elif geom_predicted.area >10000: #evaluate only the outer ring
        # pass
        # evaluate only the outer ring? # TODO issue 102
        else:
            properties[EVALUATION_FIELD_NAME] = Evaluation.TO_CHECK_NO_PREDICTION
        return properties

    @staticmethod
    def _add_multi_polygons_from_geom_to_array(geom: BaseGeometry, array):
        """
        Append valid polygons and multipolygons extracted from a given geometry to an
        existing array.

        Args:
            geom (BaseGeometry): The input geometry to process.
            array (list): An existing list to store valid polygons and multipolygons.

        Returns:
            list: A list containing valid polygons and multipolygons extracted from the
                input geometry.
        """
        if geom.is_empty or geom is None:
            # If the input geometry is empty or None, do nothing.
            pass
        else:
            # Create a GeometryCollection from the input geometry.
            geometry_collection = GeometryCollection(geom)  # noqa
            for g in geometry_collection.geoms:
                # Ensure each sub-geometry is valid.
                g = make_valid(g)
                if str(g.geom_type) in ["Polygon", "MultiPolygon"]:
                    # Append valid polygons and multipolygons to the array.
                    array.append(g)
        return array


def _calculate_geom_by_intersection_and_reference(
    geom_intersection: BaseGeometry,
    geom_reference: BaseGeometry,
    input_geometry_inner: BaseGeometry,
    is_open_domain,
    buffer_distance,
    threshold_overlap_percentage,
    threshold_exclusion_percentage,
    threshold_exclusion_area,
    threshold_inclusion_percentage,
    mitre_limit,
    partial_snapping,
    partial_snap_strategy,
    partial_snap_max_segment_length,
):
    """
    Calculates the geometry based on intersection and reference geometries.

    Args:
        geom_intersection (BaseGeometry): The intersection geometry.
        geom_reference (BaseGeometry): The reference geometry.
        is_open_domain (bool): A flag indicating whether it's a public domain
            (area not covered with reference polygon).
        threshold_exclusion_percentage (int): The threshold exclusion percentage.
        threshold_exclusion_area (int): The threshold exclusion area.
        buffer_distance (float): The buffer distance.
        threshold_overlap_percentage (int): The threshold overlap percentage.

    Returns:
        tuple: A tuple containing the resulting geometries:

        *   geom: BaseGeometry or None: The resulting geometry or None if conditions
            are not met.
        *   geom_relevant_intersection: BaseGeometry or None: The relevant
            intersection.
        *   geom_relevant_difference: BaseGeometry or None: The relevant difference.

    Notes:
        -   If the reference geometry area is 0, the overlap is set to 100%.
        -   If the overlap is less than relevant_OVERLAP_PERCENTAGE or the
            intersection area is less than relevant_OVERLAP_AREA, None is returned.
        -   Otherwise, the relevant intersection and difference geometries are
            calculated.
        -   If both relevant intersection and difference are non-empty, the final
            geometry is obtained by applying safe intersection and buffering.
        -   If only relevant intersection is non-empty, the result is the reference
            geometry.
        -   If only relevant difference is non-empty, the result is None.
    """
    od_overlap = 111  # define a specific value for defining overlap of OD
    if geom_reference.area == 0:
        overlap = od_overlap  # Open Domain

    else:
        overlap = geom_intersection.area * 100 / geom_reference.area

    if (
        overlap < threshold_exclusion_percentage
        or geom_intersection.area < threshold_exclusion_area
    ):
        return Polygon(), Polygon(), Polygon()

    if overlap >= threshold_inclusion_percentage and not overlap == od_overlap:
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
            geom_intersection, buffer_pos(geom_intersection_inner, 2 * buffer_distance)
        )

        geom_x = snap_geometry_to_reference(
            geom_x,
            geom_reference,
            max_segment_length=partial_snap_max_segment_length,
            snap_strategy=partial_snap_strategy,
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

        geom_intersection_buffered = buffer_pos(geom_intersection, 2 * buffer_distance)
        geom_difference_2 = safe_difference(geom_reference, geom_intersection_buffered)
        geom_difference_2_buffered = buffer_pos(geom_difference_2, 2 * buffer_distance)

        geom_x = safe_difference(geom_x, geom_difference_2_buffered)

        geom_x = safe_intersection(geom_x, geom_reference)

        if partial_snapping:
            geom_x = snap_geometry_to_reference(
                geom_x,
                geom_reference,
                max_segment_length=partial_snap_max_segment_length,
                snap_strategy=partial_snap_strategy,
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
            geom = _get_relevant_polygons_from_geom(geom, buffer_distance, mitre_limit)
    elif not geom_relevant_intersection.is_empty and geom_relevant_difference.is_empty:
        geom = geom_reference
    elif geom_relevant_intersection.is_empty and not geom_relevant_difference.is_empty:
        geom = geom_relevant_intersection  # (=empty geometry)
    else:
        # No relevant intersection and no relevant difference
        if is_open_domain:
            # geom = geom_relevant_intersection  # (=empty geometry)
            geom = snap_geometry_to_reference(
                geom_intersection,
                geom_reference,
                snap_strategy=partial_snap_strategy,
                tolerance=2 * buffer_distance,
                max_segment_length=partial_snap_max_segment_length,
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
        elif threshold_overlap_percentage < 0:
            # if we take a value of -1, the original border will be used
            geom = geom_intersection
        elif overlap > threshold_overlap_percentage:
            geom = geom_reference
        else:
            geom = geom_relevant_intersection  # (=empty geometry)
    return geom, geom_relevant_intersection, geom_relevant_difference


def _get_relevant_polygons_from_geom(
    geometry: BaseGeometry, buffer_distance: float, mitre_limit
):
    """
    Get only the relevant parts (polygon) from a geometry.
    Points, Lines and Polygons smaller than relevant distance are excluded from the
    result
    """
    if not geometry or geometry.is_empty:
        # If the input geometry is empty or None, do nothing.
        return geometry
    else:
        geometry = safe_unary_union(geometry)
        # Create a GeometryCollection from the input geometry.
        geometry_collection = GeometryCollection(geometry)
        array = []
        for g in geometry_collection.geoms:
            # Ensure each sub-geometry is valid.
            g = make_valid(g)
            if str(g.geom_type) in ["Polygon", "MultiPolygon"]:
                relevant_geom = buffer_neg(
                    g,
                    buffer_distance,
                    mitre_limit=mitre_limit,
                )
                if relevant_geom is not None and not relevant_geom.is_empty:
                    array.append(g)
    return safe_unary_union(array)


def _equal_geom_in_array(geom, geom_array, correction_distance, mitre_limit):
    """
    Check if a predicted geometry is equal to other predicted geometries in a list.
    Equality is defined as there is the symmetrical difference is smaller than the CORRECTION DISTANCE
    Returns True if one of the elements is equal, otherwise False
    """
    for g in geom_array:
        if geometric_equality(geom, g, correction_distance, mitre_limit):
            return True
    return False


def _calculate_inner_outer(input_geometry, relevant_distance):
    """
    calculate the inner and outer of a polygon for performance gain when using brdr_algorithm
    :param input_geometry:
    :param relevant_distance:
    :return:
    """
    input_geometry = safe_unary_union(get_parts(input_geometry))
    input_geometry_inner = buffer_neg(
        input_geometry, relevant_distance
    )  # inner part of the input that must be always available
    input_geometry_double_inner = buffer_neg(
        input_geometry, 2 * relevant_distance + MAX_OUTER_BUFFER
    )  # inner part of the input that must be always available
    # do the calculation only for the outer border of the geometry. The inner part is added afterwards
    input_geometry_outer = safe_difference(input_geometry, input_geometry_double_inner)
    return input_geometry_inner, input_geometry_outer


def _unary_union_result_dict(result_dict):
    # make a unary union for each key value in the result dict
    for key in ProcessResult.__annotations__:
        geom = result_dict.get(key, GeometryCollection())  # noqa
        if isinstance(geom, BaseGeometry) and not geom.is_empty:
            geom = safe_unary_union(geom)
        result_dict[key] = geom  # noqa
    return result_dict
