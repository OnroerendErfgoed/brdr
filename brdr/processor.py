import heapq
from abc import ABC
from abc import abstractmethod

from shapely import GeometryCollection
from shapely import LinearRing
from shapely import MultiPoint
from shapely import MultiPolygon
from shapely import Point
from shapely import Polygon
from shapely import get_parts
from shapely import make_valid
from shapely import remove_repeated_points
from shapely import segmentize
from shapely import shortest_line
from shapely.errors import GeometryTypeError
from shapely.geometry.base import BaseGeometry
from shapely.geometry.linestring import LineString
from shapely.ops import nearest_points
from shapely.ops import split

from brdr.configs import ProcessorConfig
from brdr.constants import MAX_OUTER_BUFFER
from brdr.constants import RELEVANT_DISTANCE_DECIMALS
from brdr.constants import REMARK_FIELD_NAME
from brdr.enums import OpenDomainStrategy
from brdr.enums import SnapStrategy
from brdr.geometry_utils import buffer_neg
from brdr.geometry_utils import buffer_neg_pos
from brdr.geometry_utils import buffer_pos
from brdr.geometry_utils import fill_and_remove_gaps
from brdr.geometry_utils import find_best_path_in_network
from brdr.geometry_utils import geometric_equality
from brdr.geometry_utils import get_coords_from_geometry
from brdr.geometry_utils import get_shape_index
from brdr.geometry_utils import prepare_network
from brdr.geometry_utils import safe_difference
from brdr.geometry_utils import safe_intersection
from brdr.geometry_utils import safe_symmetric_difference
from brdr.geometry_utils import safe_unary_union
from brdr.geometry_utils import safe_union
from brdr.geometry_utils import shortest_connections_between_geometries
from brdr.geometry_utils import snap_geometry_to_reference
from brdr.geometry_utils import to_multi
from brdr.logger import Logger
from brdr.typings import ProcessResult
from brdr.utils import get_relevant_polygons_from_geom
from brdr.utils import unary_union_result_dict


class BaseProcessor(ABC):
    def __init__(self, feedback, config: ProcessorConfig):
        self.logger = Logger(feedback)
        self.config = config

    @abstractmethod
    def process(self, *args, **kwargs) -> ProcessResult:
        pass

    def _postprocess_preresult(
        self,
        geom_preresult,
        geom_thematic,
        relevant_intersection,
        relevant_diff,
        relevant_distance,
        reference_union,
        mitre_limit,
        correction_distance,
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
                geom_preresult,
                buffer_pos(geom_thematic, correction_distance),
            )
            result_diff_min = safe_difference(
                geom_thematic,
                buffer_pos(geom_preresult, correction_distance),
            )
            result_diff = safe_unary_union([result_diff_plus, result_diff_min])
            return unary_union_result_dict(
                {
                    "result": geom_preresult,
                    "result_diff": result_diff,
                    "result_diff_plus": result_diff_plus,
                    "result_diff_min": result_diff_min,
                    "result_relevant_intersection": relevant_intersection,
                    "result_relevant_diff": relevant_diff,
                    "properties": {REMARK_FIELD_NAME: remark},
                }
            )
        # Process array

        buffer_distance = relevant_distance / 2
        result = []
        geom_thematic_for_add_delete = geom_thematic

        if self.config.od_strategy == OpenDomainStrategy.EXCLUDE:
            geom_thematic_for_add_delete = safe_intersection(
                geom_thematic_for_add_delete, reference_union
            )
            geom_preresult = safe_intersection(geom_preresult, reference_union)

        if not (geom_thematic is None or geom_thematic.is_empty):
            # Correction for circles
            # calculate ratio to see if it is a circle, and keep the original geometry
            #  if a circle: (Polsby-Popper score)
            if (
                get_shape_index(geom_thematic.area, geom_thematic.length)
                > self.config.threshold_circle_ratio
            ):
                remark = "Circle detected: -->resulting geometry = original geometry"
                self.logger.feedback_debug(remark)
                return unary_union_result_dict(
                    {"result": geom_thematic, "properties": {REMARK_FIELD_NAME: remark}}
                )

            # Correction for unchanged geometries
            # if safe_symmetric_difference(geom_preresult, geom_thematic).is_empty:
            if geometric_equality(
                geom_preresult,
                geom_thematic,
                correction_distance=correction_distance,
                mitre_limit=mitre_limit,
            ):
                remark = "Unchanged geometry: -->resulting geometry = original geometry"
                self.logger.feedback_debug(remark)
                return unary_union_result_dict(
                    {"result": geom_thematic, "properties": {REMARK_FIELD_NAME: remark}}
                )

        # Corrections for areas that differ more than the relevant distance
        geom_thematic_dissolved = buffer_pos(
            buffer_neg(
                buffer_pos(
                    geom_preresult,
                    correction_distance,
                    mitre_limit=mitre_limit,
                ),
                2 * correction_distance,
                mitre_limit=mitre_limit,
            ),
            correction_distance,
            mitre_limit=mitre_limit,
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
                    mitre_limit=mitre_limit,
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
                    mitre_limit=mitre_limit,
                ),
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
                safe_symmetric_difference(geom_thematic_result, geom_thematic),
                correction_distance,
                mitre_limit=mitre_limit,
            ),
            correction_distance,
            mitre_limit=mitre_limit,
        )
        geom_result_diff_plus = buffer_pos(
            buffer_neg(
                safe_difference(
                    geom_thematic_result,
                    geom_thematic,
                ),
                correction_distance,
                mitre_limit=mitre_limit,
            ),
            correction_distance,
            mitre_limit=mitre_limit,
        )
        geom_result_diff_min = buffer_pos(
            buffer_neg(
                safe_difference(
                    geom_thematic,
                    geom_thematic_result,
                ),
                correction_distance,
                mitre_limit=mitre_limit,
            ),
            correction_distance,
            mitre_limit=mitre_limit,
        )

        return unary_union_result_dict(
            {
                "result": geom_thematic_result,
                "result_diff": geom_result_diff,
                "result_diff_plus": geom_result_diff_plus,
                "result_diff_min": geom_result_diff_min,
                "result_relevant_intersection": relevant_intersection,
                "result_relevant_diff": relevant_diff,
                "properties": {REMARK_FIELD_NAME: remark},
            }
        )

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
            > relevant_distance * self.config.partial_snap_max_segment_length * 2
            # * factor
            # * factor
            # * 4  # There could be a better way to exclude invalid connection-lines?
        ):
            return LineString()

        return connection_line

    def _create_virtual_reference(
        self,
        geometry,
        relevant_distance,
        reference_union,
        correction_distance,
        mitre_limit,
        use_outer_boundary=False,
    ):
        """
        Functions that creates a 'virtual reference polygon' for areas that are not covered by reference-polygons.
        :param geometry:
        :param relevant_distance:
        :param use_outer_boundary: when outer is True, the outer boundary is used, inner is not used
        :return:
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
        )  # Both OD-parts are SNAPPED added
        if (
            use_outer_boundary
        ):  # when outer is True, the outer boundary is used, inner is not used
            geom_1 = safe_difference(geometry, virtual_reference)
            geom_2 = buffer_neg_pos(geom_1, buffer_distance)
            geom_3 = safe_intersection(geom_2, geometry)
            virtual_reference = safe_unary_union([geom_3, virtual_reference])
        return virtual_reference

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


class SnapGeometryProcessor(BaseProcessor):
    def process(
        self,
        *,
        correction_distance,
        input_geometry: BaseGeometry,
        mitre_limit,
        ref_intersections_geoms,
        reference_union,
        relevant_distance,
        snap_max_segment_length,
        snap_strategy,
    ) -> ProcessResult:
        snapped = []
        virtual_reference = Polygon()
        if self.config.od_strategy != OpenDomainStrategy.EXCLUDE:
            virtual_reference = self._create_virtual_reference(
                input_geometry,
                relevant_distance,
                reference_union,
                correction_distance,
                mitre_limit,
                # TODO missing params
                False,
            )
        if self.config.od_strategy == OpenDomainStrategy.EXCLUDE:
            pass
        elif self.config.od_strategy == OpenDomainStrategy.AS_IS:
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
        geom_preresult = safe_unary_union(snapped)
        result_dict = self._postprocess_preresult(
            geom_preresult,
            input_geometry,
            GeometryCollection(),
            GeometryCollection(),
            relevant_distance,
            reference_union,
            mitre_limit,
            correction_distance,
        )
        return result_dict


class DieussaertGeometryProcessor(BaseProcessor):
    def process(
        self,
        input_geometry,
        input_geometry_outer,
        input_geometry_inner,
        ref_intersections_geoms,
        relevant_distance,
        mitre_limit,
        reference_union,
        correction_distance,
    ) -> ProcessResult:
        buffer_distance = relevant_distance / 2
        (
            preresult,
            relevant_intersection_array,
            relevant_diff_array,
        ) = self._calculate_intersection_between_geometry_and_od(
            input_geometry_outer,
            input_geometry_inner,
            relevant_distance,
            reference_union,
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
                geom_intersection,
                geom_reference,
                input_geometry_inner,
                False,
                buffer_distance,
                mitre_limit,
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
        result_dict = self._postprocess_preresult(
            geom_preresult,
            input_geometry,
            relevant_intersection,
            relevant_diff,
            relevant_distance,
            reference_union,
            mitre_limit,
            correction_distance,
        )
        return result_dict

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
                geom_intersection,
                geom_reference,
                input_geometry_inner,
                True,
                buffer_distance,
                mitre_limit,
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
        geom_intersection: BaseGeometry,
        geom_reference: BaseGeometry,
        input_geometry_inner: BaseGeometry,
        is_open_domain,
        buffer_distance,
        mitre_limit,
    ):
        """
        Calculates the geometry based on intersection and reference geometries.

        Args:
            geom_intersection (BaseGeometry): The intersection geometry.
            geom_reference (BaseGeometry): The reference geometry.
            is_open_domain (bool): A flag indicating whether it's a public domain
                (area not covered with reference polygon).
            buffer_distance (float): The buffer distance.

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
    def process(
        self,
        *,
        input_geometry,
        reference_elements,
        reference_union,
        mitre_limit,
        correction_distance,
        relevant_distance=1,
        snap_strategy: SnapStrategy = SnapStrategy.NO_PREFERENCE,
    ) -> ProcessResult:
        input_geometry = to_multi(input_geometry)
        input_geometry_buffered = buffer_pos(
            input_geometry,
            relevant_distance * self.config.buffer_multiplication_factor,
        )
        reference = safe_unary_union(
            safe_intersection(reference_elements, input_geometry_buffered)
        )
        geom_processed_list = []
        if isinstance(input_geometry, MultiPolygon):
            for polygon in input_geometry.geoms:
                exterior = polygon.exterior
                interiors = polygon.interiors
                exterior_processed = self._process_by_network(
                    exterior,
                    reference,
                    relevant_distance,
                    close_output=True,
                )
                interiors_processed = []
                for i in interiors:
                    i_processed = self._process_by_network(
                        i,
                        reference,
                        relevant_distance,
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
            reference_union,
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
                input_geometry = geom_to_process,
                reference_intersection = reference_intersection,
                reference_coords_intersection = reference_coords_intersection,
                thematic_difference = thematic_difference,
                relevant_distance = relevant_distance,
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
        geom_to_process_segmentized = segmentize(
            geom_to_process_line, self.config.partial_snap_max_segment_length
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

    def _get_processed_network_path(
        self,
        input_geometry,
        reference_intersection,
        reference_coords_intersection,
        thematic_difference,
        relevant_distance,
    ):
        thematic_points = self._get_thematic_points(
            input_geometry, reference_intersection
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
        )  # removed as we first will filter these lines

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
            input_geometry, network, self.config.snap_strategy, relevant_distance
        )
        # if geom_processed is None:
        #     # add original so a connected path will be found
        #     segments.append(geom_to_process)
        #     network = prepare_network(segments)
        #     geom_processed = find_longest_path_in_network(
        #         geom_to_process, network, snap_strategy, relevant_distance
        #     )

        return geom_processed


class AlignerGeometryProcessor(BaseProcessor):
    def process(
        self,
        mitre_limit,
        reference_union,
        reference_elements,
        correction_distance,
        reference_items,
        reference_tree,
        dict_reference,
        input_geometry: BaseGeometry,
        relevant_distance: int = 1,
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

        self.config.od_strategy = od_strategy
        self.config.threshold_overlap_percentage = threshold_overlap_percentage

        if isinstance(input_geometry, GeometryCollection):
            raise ValueError(
                "GeometryCollection as input is not supported. Please use the individual geometries from the GeometryCollection as input."
            )
        elif isinstance(input_geometry, (Polygon, MultiPolygon)):
            # Processing thematic polygons
            if self.config.area_limit and input_geometry.area > self.config.area_limit:
                message = f"The input polygon is too large to process: input area {str(input_geometry.area)} m², limit area: {str(self.config.area_limit)} m²."
                raise ValueError(message)
            self.logger.feedback_debug("process geometry")

            # For calculations with RD=0 the original input is returned
            if relevant_distance == 0:
                self.logger.feedback_debug("Calculation for RD = 0")
                return unary_union_result_dict(
                    {
                        "result": input_geometry,
                        "properties": {
                            REMARK_FIELD_NAME: "relevant distance 0 --> original geometry returned"
                        },
                    }
                )

            # CALCULATE INNER and OUTER INPUT GEOMETRY for performance optimisation on big geometries
            # combine all parts of the input geometry to one polygon
            input_geometry_inner, input_geometry_outer = self._calculate_inner_outer(
                input_geometry, relevant_distance
            )
            # get a list of all ref_ids that are intersecting the thematic geometry; we take it bigger because we want to check if there are also reference geometries on a relevant distance.
            input_geometry_outer_buffered = buffer_pos(
                input_geometry_outer,
                relevant_distance * self.config.buffer_multiplication_factor,
            )
            ref_intersections = reference_items.take(
                reference_tree.query(input_geometry_outer_buffered)
            ).tolist()

            ref_intersections_geoms = []
            all_polygons = True
            for key_ref in ref_intersections:
                ref_geom = dict_reference[key_ref]
                ref_intersections_geoms.append(ref_geom)
                if not isinstance(ref_geom, (Polygon, MultiPolygon)):
                    all_polygons = False

            if all_polygons:
                processor = DieussaertGeometryProcessor(
                    self.logger.feedback, self.config
                )
                return processor.process(
                    input_geometry,
                    input_geometry_outer,
                    input_geometry_inner,
                    ref_intersections_geoms,
                    relevant_distance,
                    mitre_limit,
                    reference_union,
                    correction_distance,
                )
        processor = NetworkGeometryProcessor(self.logger.feedback, self.config)
        return processor.process(
            input_geometry=input_geometry,
            reference_elements=reference_elements,
            reference_union=reference_union,
            relevant_distance=relevant_distance,
            snap_strategy=self.config.partial_snap_strategy,
            mitre_limit=mitre_limit,
            correction_distance=correction_distance,
        )

    @staticmethod
    def _calculate_inner_outer(input_geometry, relevant_distance):
        """
        calculate the inner and outer of a polygon for performance gain when using dieussaert_algorithm
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
        # do the calculation only for the outer border of the geometry. The inner part is added afterward
        input_geometry_outer = safe_difference(
            input_geometry, input_geometry_double_inner
        )
        return input_geometry_inner, input_geometry_outer
