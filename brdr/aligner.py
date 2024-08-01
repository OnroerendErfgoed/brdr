import json
import logging
import os
from collections import defaultdict
from math import pi
from typing import Iterable

import numpy as np
from shapely import GeometryCollection
from shapely import Polygon
from shapely import STRtree
from shapely import equals
from shapely import get_exterior_ring
from shapely import get_interior_ring
from shapely import get_num_interior_rings
from shapely import get_parts
from shapely import make_valid
from shapely import polygons
from shapely import remove_repeated_points
from shapely import to_geojson
from shapely import unary_union
from shapely.geometry.base import BaseGeometry

from brdr.constants import BUFFER_MULTIPLICATION_FACTOR
from brdr.constants import CORR_DISTANCE
from brdr.constants import DEFAULT_CRS
from brdr.constants import DOWNLOAD_LIMIT
from brdr.constants import THRESHOLD_CIRCLE_RATIO
from brdr.enums import GRBType
from brdr.enums import OpenbaarDomeinStrategy
from brdr.geometry_utils import buffer_neg
from brdr.geometry_utils import buffer_neg_pos
from brdr.geometry_utils import buffer_pos
from brdr.geometry_utils import calculate_geom_by_intersection_and_reference
from brdr.geometry_utils import safe_difference
from brdr.geometry_utils import safe_intersection
from brdr.geometry_utils import safe_symmetric_difference
from brdr.geometry_utils import safe_union
from brdr.loader import DictLoader
from brdr.loader import GRBActualLoader
from brdr.loader import GeoJsonFileLoader
from brdr.loader import GeoJsonLoader
from brdr.loader import GeoJsonUrlLoader
from brdr.loader import Loader
from brdr.logger import Logger
from brdr.typings import ProcessResult
from brdr.utils import diffs_from_dict_series, dict_predicted_by_keys
from brdr.utils import geojson_from_dict
from brdr.utils import get_series_geojson_dict
from brdr.utils import get_breakpoints_zerostreak
from brdr.utils import get_collection
from brdr.utils import merge_process_results
from brdr.utils import write_geojson

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(message)s", datefmt="%d-%b-%y %H:%M:%S"
)


###################


class Aligner:
    """
    This class is used to compare the thematic data with the reference data.
    The reference data can be loaded in different ways, for example by using the GRB data.
    The thematic data can be loaded by using a geojson file.
    The class can be used to compare the thematic data with the reference data.
    """

    def __init__(
        self,
        *,
        feedback=None,
        relevant_distance=1,
        threshold_overlap_percentage=50,
        od_strategy=OpenbaarDomeinStrategy.SNAP_SINGLE_SIDE,
        crs=DEFAULT_CRS,
    ):
        """
        Initializes the Aligner object

        Args:
            feedback (object, optional): Feedback object that can be added to show
                feedback in QGIS. Defaults to None.
            relevant_distance (int, optional): The relevant distance (in meters) for
                processing. Defaults to 1.
            od_strategy (int, optional): The strategy to determine how to handle
                information outside the reference polygons (Openbaar Domein)
                (default 1: SNAP_SINGLE_SIDE)
            threshold_overlap_percentage (int, optional): Threshold (%) to determine
                from which overlapping-percentage a reference-polygon has to be included
                when there aren't relevant intersections or relevant differences
                (default 50%)

        """
        self.logger = Logger(feedback)
        self.relevant_distance = relevant_distance
        self.od_strategy = od_strategy
        self.threshold_overlap_percentage = threshold_overlap_percentage

        # PROCESSING DEFAULTS
        # thematic
        # name of the identifier-field of the thematic data (id has to be unique)
        self.name_thematic_id = "theme_identifier"
        # dictionary to store all thematic geometries to handle
        self.dict_thematic: dict[str, BaseGeometry] = {}

        # reference
        self.name_reference_id = "ref_identifier"  # name of the identifier-field of the reference data (id has to be unique,f.e CAPAKEY for GRB-parcels)
        self.dict_reference = {}  # dictionary to store all reference geometries
        self.reference_union = None  # to save a unioned geometry of all reference polygons; needed for calculation in most OD-strategies

        # output-dictionaries (when processing dict_thematic)
        self.dict_result: dict[str, ProcessResult] = {}
        # dictionary with the 'predicted' results, grouped by relevant distance
        self.dict_predicted = dict[float, dict[str, ProcessResult]]

        # Coordinate reference system
        # thematic geometries and reference geometries are assumed to be in the same CRS
        # before loading into the Aligner. No CRS-transformation will be performed.
        # When loading data, CRS is expected to be a projected CRS with units in 'meter (m)'.
        # Default EPSG:31370 (Lambert72), alternative: EPSG:3812 (Lambert2008)
        self.CRS = crs

        self.logger.feedback_info("Aligner initialized")

    def buffer_distance(self):
        return self.relevant_distance / 2

    def process_geometry(
        self,
        geometry: BaseGeometry,
        relevant_distance=1,
        od_strategy=OpenbaarDomeinStrategy.SNAP_SINGLE_SIDE,
        threshold_overlap_percentage=50,
    ) -> ProcessResult:
        """
        method to align a geometry to the reference layer

        Args:
            geometry (BaseGeometry): The input geometric object.
            relevant_distance
            od_strategy
            threshold_overlap_percentage (int): The buffer distance (positive or negative).

        Returns:
            ProcessResult : A dict containing the resulting geometries:

            *   result (BaseGeometry): The resulting output geometry
            *   result_diff (BaseGeometry): The resulting difference output geometry
            *   result_diff_plus (BaseGeometry): The resulting positive difference output
                geometry
            *   result_diff_min (BaseGeometry): The resulting negative difference output
                geometry
            *   relevant_intersection (BaseGeometry): The relevant_intersection
            *   relevant_difference (BaseGeometry): The relevant_difference
        Notes:
            -
        Example:
        """
        self.logger.feedback_debug("process geometry")
        self.relevant_distance = relevant_distance
        self.od_strategy = od_strategy
        self.threshold_overlap_percentage = threshold_overlap_percentage
        # array with all relevant parts of a thematic geometry; initial empty Polygon
        preresult = [Polygon()]
        (
            geometry,
            preresult,
            relevant_intersection_array,
            relevant_diff_array,
        ) = self._calculate_intersection_between_geometry_and_od(geometry, preresult)
        # get a list of all ref_ids that are intersecting the thematic geometry
        ref_intersections = self.reference_items.take(
            self.reference_tree.query(geometry)
        ).tolist()
        for key_ref in ref_intersections:
            geom_reference = self.dict_reference[key_ref]
            geom_intersection = safe_intersection(geometry, geom_reference)
            if geom_intersection.is_empty or geom_intersection is None:
                continue
            self.logger.feedback_debug("calculate intersection")
            (
                geom,
                relevant_intersection,
                relevant_diff,
            ) = calculate_geom_by_intersection_and_reference(
                geom_intersection,
                geom_reference,
                False,
                self.relevant_distance / 2,
                self.threshold_overlap_percentage,
            )
            self.logger.feedback_debug("intersection calculated")
            preresult = self.add_multi_polygons_from_geom_to_array(geom, preresult)
            relevant_intersection_array = self.add_multi_polygons_from_geom_to_array(
                relevant_intersection, relevant_intersection_array
            )
            relevant_diff_array = self.add_multi_polygons_from_geom_to_array(
                relevant_diff, relevant_diff_array
            )
        # UNION INTERMEDIATE LAYERS
        relevant_intersection = unary_union(relevant_intersection_array)
        if relevant_intersection is None or relevant_intersection.is_empty:
            relevant_intersection = Polygon()
        relevant_diff = unary_union(relevant_diff_array)
        if relevant_diff is None or relevant_diff.is_empty:
            relevant_diff = Polygon()

        # POSTPROCESSING
        result_dict = self._postprocess_preresult(preresult, geometry)

        result_dict["result_relevant_intersection"] = relevant_intersection
        result_dict["result_relevant_diff"] = relevant_diff
        return result_dict

    def process_dict_thematic(
        self,
        relevant_distance=1,
        od_strategy=OpenbaarDomeinStrategy.SNAP_SINGLE_SIDE,
        threshold_overlap_percentage=50,
    ) -> dict[str, ProcessResult]:
        """
        Aligns a thematic dictionary of geometries to the reference layer based on
        specified parameters. - method to align a thematic dictionary to the reference
        layer

        Args:
            relevant_distance (float, optional): The relevant distance (in meters) for
                processing. Defaults to 1.
            od_strategy (int, optional): The strategy for overlap detection.
                Defaults to 1.
            threshold_overlap_percentage (float, optional): The threshold percentage for
                considering full overlap. Defaults to 50.

        Returns:
            dict: A dict containing processed data for each thematic key:
                - result: Aligned thematic data.
                - result_diff: global differences between thematic data and reference data.
                - result_diff_plus: Positive differences.
                - result_diff_min: Negative differences.
                - relevant_intersection: relevant intersections.
                - relevant_diff: relevant differences.

        """
        self.dict_result = {}
        for key in self.dict_thematic:
            self.logger.feedback_info("thematic id to process: " + str(key))
            self.dict_result[key] = self.process_geometry(
                self.dict_thematic[key],
                relevant_distance,
                od_strategy,
                threshold_overlap_percentage,
            )
        return self.dict_result

    def predictor(
        self,
        relevant_distances=np.arange(0, 300, 10, dtype=int) / 100,
        od_strategy=OpenbaarDomeinStrategy.SNAP_SINGLE_SIDE,
        threshold_overlap_percentage=50,
    ):
        """
        Predicts the 'most interesting' relevant distances for changes in thematic elements based on a distance series.

        This function analyzes a set of thematic geometries (`self.dict_thematic`) to identify potentially
        interesting distances where changes occur. It performs the following steps:

        1. **Process Distance Series:**
            - Calculates a series of results for different distances specified by `relevant_distances`.
            - This calculation might involve functions like `self.process_series` (implementation details likely depend on your specific code).

        2. **Calculate Difference Metrics:**
            - Analyzes the results from the distance series to compute difference metrics
              between thematic elements at each distance (using `diffs_from_dict_series`).

        3. **Identify Breakpoints and Zero-Streaks:**
            - For each thematic geometry, it identifies potential "breakpoints" where the difference metric changes sign (from positive to negative or vice versa).
            - It also identifies "zero-streaks" which are consecutive distances with a difference metric close to zero (potentially indicating minimal change).

        4. **Predict Interesting Distances:**
            - The function considers distances corresponding to breakpoints and zero-streaks as potentially interesting for further analysis.
            - These distances are stored in a dictionary (`dict_predicted`) with the thematic element key as the outer key.
            - Additionally, the corresponding results from the distance series for those distances are included.

        5. **Filter Results:**
            - The function might further filter the predicted results for each thematic element based on the element key (using `filter_resulting_series_by_key`).

        Args:
            relevant_distances (np.ndarray, optional): A NumPy array of distances to be analyzed. Defaults to np.arange(0.1, 5.05, 0.1).
            od_strategy (OpenbaarDomeinStrategy, optional): A strategy for handling open data in the processing (implementation specific). Defaults to OpenbaarDomeinStrategy.SNAP_SINGLE_SIDE.
            threshold_overlap_percentage (int, optional): A percentage threshold for considering full overlap in the processing (implementation specific). Defaults to 50.

        Returns:
            dict: A dictionary containing predicted interesting distances for each thematic element.
                - Keys: Thematic element identifiers from `self.dict_thematic`.
                - Values: Dictionaries with the following structure for each thematic element:
                    - Keys: Distances identified as interesting (breakpoints or zero-streaks).
                    - Values: dicts containing results (likely specific to your implementation) from the distance series for the corresponding distance.

        Logs:
            - Debug logs the thematic element key being processed.
        """
        dict_predicted = defaultdict(dict)
        dict_series = self.process_series(
            relevant_distances=relevant_distances,
            od_strategy=od_strategy,
            threshold_overlap_percentage=threshold_overlap_percentage,
        )

        diffs_dict = diffs_from_dict_series(dict_series, self.dict_thematic)

        for theme_id, diffs in diffs_dict.items():
            if len(diffs) != len(relevant_distances):
                logging.warning(
                    f"Number of computed diffs for thematic element {theme_id} does not match the number of relevant distances."
                )
                continue
            diff_values = list(diffs.values())
            breakpoints, zero_streaks = get_breakpoints_zerostreak(
                relevant_distances, diff_values
            )
            logging.debug(str(theme_id))
            for zs in zero_streaks:
                dict_predicted[zs[0]][theme_id] = dict_series[zs[0]][theme_id]

        self.dict_predicted = dict_predicted


        return dict_series, dict_predicted, diffs_dict,

    def process_series(
        self,
        relevant_distances: Iterable[float],
        od_strategy=OpenbaarDomeinStrategy.SNAP_SINGLE_SIDE,
        threshold_overlap_percentage=50,
    ) -> dict[float, dict[str, ProcessResult]]:
        """
        Calculates the resulting dictionaries for thematic data based on a series of relevant
            distances.

        Args:
            relevant_distances (Iterable[float]): A series of relevant distances (in meters) to
                process
            od_strategy (int, optional): The strategy for overlap detection.
                Defaults to 1.
            threshold_overlap_percentage (float, optional): The threshold percentage for
                considering full overlap. Defaults to 50.

        Returns:
            dict: A dictionary containing the resulting dictionaries for a series of relevant distances:

                {
                    'relevant_distance_1': {theme_id_1: (ProcessResult), theme_id_2: (ProcessResult), ...},
                    'relevant_distance_2': {theme_id_1: (ProcessResult), theme_id_2: (ProcessResult), ...},
                    ...
                }
        """
        self.logger.feedback_debug("Process series" + str(relevant_distances))
        self.od_strategy = od_strategy
        self.threshold_overlap_percentage = threshold_overlap_percentage
        dict_series = {}
        for s in relevant_distances:
            self.logger.feedback_info(
                "Processing series - relevant_distance (m):"
                + str(s)
                + " with ODStrategy "
                + str(self.od_strategy)
            )
            dict_series[s] = self.process_dict_thematic(s, od_strategy)
        self.logger.feedback_info(
            "End of processing series: " + str(relevant_distances)
        )
        return dict_series

    def get_formula(self, geometry: BaseGeometry, with_geom=False):
        """
        Calculates formula-related information based on the input geometry.

        Args:
            geometry (shapely.geometry object): The input geometry.
            with_geom (bool, optional): Whether to include geometry information in the
                output. Defaults to False.

        Returns:
            dict: A dictionary containing formula-related data:

                -   'full': True if the intersection is the same as the reference
                    geometry, else False.
                -   'area': Area of the intersection or reference geometry.
                -   'percentage': Percentage of intersection area relative to the
                    reference geometry.
                -   'geometry': GeoJSON representation of the intersection (if
                    with_geom is True).
        """
        dict_formula = {}
        ref_intersections = self.reference_items.take(
            self.reference_tree.query(geometry)
        ).tolist()
        for key_ref in ref_intersections:
            geom_reference = self.dict_reference[key_ref]
            geom_intersection = make_valid(safe_intersection(geometry, geom_reference))
            if geom_intersection.is_empty or geom_intersection is None:
                continue
            if equals(geom_intersection, geom_reference):
                full = True
                area = round(geom_reference.area, 2)
                perc = 100
                geom = to_geojson(geom_reference)
            else:
                perc = round(geom_intersection.area * 100 / geom_reference.area, 2)
                if perc < 0.01:
                    continue
                elif perc > 99.99:
                    full = True
                    area = round(geom_reference.area, 2)
                    perc = 100
                    geom = to_geojson(geom_reference)
                else:
                    full = False
                    area = round(geom_intersection.area, 2)
                    geom = to_geojson(geom_intersection)
            if not with_geom:
                geom = None

            dict_formula[key_ref] = {
                "full": full,
                "area": area,
                "percentage": perc,
                "geometry": geom,
            }

        self.logger.feedback_debug(str(dict_formula))
        return dict_formula

    def get_results_as_dict(self, merged=True):
        """
        get a dict of the results

        Args:
            merged (bool, optional): Whether to merge the results for each thematic element. Defaults to True.
        """
        if merged:
            return merge_process_results(self.dict_result)
        return self.dict_result

    def get_results_as_geojson(self, formula=False, merged=True):
        """
        convert the results to geojson feature collections

        Args:
            formula (bool, optional): Whether to include formula-related information in the output. Defaults to False.
            merged (bool, optional): Whether to merge the results for each thematic element. Defaults to True.
        """
        results_dict = self.dict_result
        if merged:
            results_dict = merge_process_results(results_dict)

        return self.get_predictions_as_geojson(
            formula,
            {self.relevant_distance: results_dict},
        )

    def get_predictions_as_dict(self):
        """
        get a dict with results for the predicted relevant distances
        """
        return self.dict_predicted

    def get_predictions_as_geojson(self, formula=False, series_dict=None):
        """
        get a dictionary containing of the resulting geometries as geojson, based on the 'predicted' relevant distances.
        Optional: The descriptive formula is added as an attribute to the result"""

        series_dict = series_dict or self.dict_predicted
        prop_dictionary = defaultdict(dict)

        for relevant_distance, results_dict in series_dict.items():
            for theme_id, process_results in results_dict.items():
                result = process_results["result"]
                if formula:
                    formula = self.get_formula(result)
                    prop_dictionary[relevant_distance][theme_id] = {
                        "formula": json.dumps(formula)
                    }

        return get_series_geojson_dict(
            series_dict,
            crs=self.CRS,
            id_field=self.name_thematic_id,
            series_prop_dict=prop_dictionary,
        )

    def get_reference_as_geojson(self):
        """
        get a geojson of the reference polygons
        """
        return geojson_from_dict(
            self.dict_reference, self.CRS, self.name_reference_id, geom_attributes=False
        )

    def export_results(self, path, formula=True):
        """
        Exports analysis results as GeoJSON files.

        This function exports 6 GeoJSON files containing the analysis results to the specified `path`.

        Args:
            path (str): The path to the directory where the GeoJSON files will be saved.
            formula (bool, optional): Whether to include formula-related information in the output. Defaults to True.

        Details of exported files:
            - result.geojson: Contains the original thematic data from `self.dict_result`.
            - result_diff.geojson: Contains the difference between the original and predicted data from `self.dict_result_diff`.
            - result_diff_plus.geojson: Contains results for areas that are added (increased area).
            - result_diff_min.geojson: Contains results for areas that are removed (decreased area).
            - result_relevant_intersection.geojson: Contains the areas with relevant intersection that has to be included in the result.
            - result_relevant_difference.geojson: Contains the areas with relevant difference that has to be excluded from the result.
        """
        fcs = self.get_results_as_geojson(formula)
        for name, fc in fcs.items():
            write_geojson(os.path.join(path, name + ".geojson"), fc)

    def _prepare_reference_data(self):
        """
        Prepares reference data for spatial queries and analysis.

        It performs the following tasks:

        1. **Optimizes spatial queries:**
            - Creates a Spatial Relationship Tree (STRtree) using `STRtree` for efficient spatial queries against the reference data in `self.dict_reference`.
            - Converts the dictionary keys (reference identifiers) to a NumPy array for potential performance benefits in future operations.

        2. **Clears reference union:**
            - Sets `self.reference_union` to `None`. This variable stores the combined geometry of all reference data,
              and it's cleared here to indicate that it needs to be recalculated if requested later.

        Returns:
            None
        """
        # create an SRTree for performance optimisation
        self.logger.feedback_info(
            "length of reference_dict: " + str(len(self.dict_reference))
        )
        self.reference_tree = STRtree(list(self.dict_reference.values()))
        self.reference_items = np.array(list(self.dict_reference.keys()))
        # clear the reference_union, so it will be recalculated on request when needed
        self.reference_union = None
        return

    def _calculate_intersection_between_geometry_and_od(self, geometry, preresult):
        # Calculate the intersection between thematic and Openbaar Domein
        relevant_intersection_array = [Polygon()]
        relevant_difference_array = [Polygon()]
        geom_thematic_od = Polygon()

        if self.od_strategy == OpenbaarDomeinStrategy.EXCLUDE:
            # Completely exclude everything that is not on the reference layer
            self.logger.feedback_debug("OD-strategy EXCLUDE")
            # Remove from the thematic layer all parts that are not on the reference layer
            # !!this strategy adapts the input-geometry!!
            geometry = safe_intersection(geometry, self._get_reference_union())
        elif self.od_strategy == OpenbaarDomeinStrategy.AS_IS:
            # All parts that are not covered by the reference layer are added to the
            #         resulting geometry AS IS
            self.logger.feedback_debug("OD-strategy AS IS")
            # all OD-parts wil be added AS IS
            geom_thematic_od = safe_difference(geometry, self._get_reference_union())
        elif self.od_strategy == OpenbaarDomeinStrategy.SNAP_SINGLE_SIDE:
            # Everything that falls within the relevant distance over the
            #  plot boundary is snapped to the plot.
            #  Only the inner-reference-boundaries are used.
            #  The outer-reference-boundaries are not used.
            self.logger.feedback_debug("OD-strategy SNAP_SINGLE_SIDE")
            # geom of OD
            geom_od = safe_difference(geometry, self._get_reference_union())
            # only the relevant parts of OD
            geom_od_neg_pos = buffer_neg_pos(geom_od, self.buffer_distance())
            # geom_thematic_od = safe_intersection(geom_od_neg_pos,geom_od)# resulting
            # thematic OD
            geom_od_neg_pos_buffered = buffer_pos(
                geom_od_neg_pos, self.buffer_distance()
            )  # include parts
            geom_thematic_od = safe_intersection(
                geom_od_neg_pos_buffered, geom_od
            )  # resulting thematic OD
        elif self.od_strategy == OpenbaarDomeinStrategy.SNAP_ALL_SIDE:
            # Everything that falls within the relevant distance over
            # the plot boundary is snapped to the plot.
            #  Inner-reference-boundaries and outer-reference-boundaries are used.
            self.logger.feedback_debug("OD-strategy SNAP BOTH SIDED")
            (
                geom_thematic_od,
                relevant_difference_array,
                relevant_intersection_array,
            ) = self._od_snap_all_side(
                geometry,
                relevant_difference_array,
                relevant_intersection_array,
            )
        elif self.od_strategy == OpenbaarDomeinStrategy.SNAP_FULL_AREA_SINGLE_SIDE:
            # Strategy useful for bigger areas.
            # integrates the entire inner area of the input geometry,
            # so Openbaar Domein of the inner area is included in the result
            # Combines SNAP_SINGLE_SIDE with the inner area
            self.logger.feedback_debug(
                "OD-strategy Full-area-variant of OD-SNAP_SINGLE_SIDE"
            )
            geom_thematic_od = self._od_full_area(geometry)
        elif self.od_strategy == OpenbaarDomeinStrategy.SNAP_FULL_AREA_ALL_SIDE:
            # Strategy useful for bigger areas.
            # integrates the entire inner area of the input geometry,
            # so Openbaar Domein of the inner area is included in the result
            # Combines SNAP_ALL_SIDE with the inner area
            self.logger.feedback_debug(
                "OD-strategy Full-area-variant of OD-SNAP_ALL_SIDE"
            )
            # first part is a copy of OD_ALL_SIDE
            (
                geom_thematic_od,
                relevant_difference_array,
                relevant_intersection_array,
            ) = self._od_snap_all_side(
                geometry,
                relevant_difference_array,
                relevant_intersection_array,
            )
            # This part is a copy of SNAP_FULL_AREA_SINGLE_SIDE
            geom_theme_od_min_clipped_plus_buffered_clipped = self._od_full_area(
                geometry
            )
            # UNION the calculation of  OD-SNAP_ALL_SIDE with FULL AREA of OD-SNAP_FULL_AREA_SINGLE_SIDE
            geom_thematic_od = safe_union(
                geom_theme_od_min_clipped_plus_buffered_clipped, geom_thematic_od
            )

        elif self.od_strategy == OpenbaarDomeinStrategy.SNAP_SINGLE_SIDE_VARIANT_1:
            # OD-strategy SNAP_SINGLE_SIDE - variant 1:
            # Everything that falls within the relevant distance over the
            #  plot boundary is snapped to the plot.
            #  Only the inner-reference-boundaries are used.
            #  The outer-reference-boundaries are not used.
            self.logger.feedback_debug("OD-strategy SNAP_SINGLE_SIDE - variant 1")
            # geom of OD
            geom_od = safe_difference(geometry, self._get_reference_union())
            # only the relevant parts of OD
            geom_od_neg_pos = buffer_neg_pos(geom_od, self.buffer_distance())
            # resulting thematic OD
            geom_thematic_od = safe_intersection(geom_od_neg_pos, geom_od)
        elif self.od_strategy == OpenbaarDomeinStrategy.SNAP_SINGLE_SIDE_VARIANT_2:
            # OD-strategy SNAP_SINGLE_SIDE - variant 2:
            # Everything that falls within the relevant distance over the
            #  plot boundary is snapped to the plot.
            #  Only the inner-reference-boundaries are used.
            #  The outer-reference-boundaries are not used.
            self.logger.feedback_debug("OD-strategy SNAP_SINGLE_SIDE - variant 2")
            # TODO: Drop?
            pass

        # ADD THEMATIC_OD
        preresult = self.add_multi_polygons_from_geom_to_array(
            geom_thematic_od, preresult
        )
        return (
            geometry,
            preresult,
            relevant_intersection_array,
            relevant_difference_array,
        )

    def _od_full_area(self, geometry):
        geom_theme_od = safe_difference(geometry, self._get_reference_union())
        geom_theme_min_buffered = buffer_neg(
            buffer_pos(
                buffer_neg(geometry, self.relevant_distance),
                self.buffer_distance(),
            ),
            self.buffer_distance(),
        )
        geom_theme_od_clipped_min_buffered = safe_intersection(
            geom_theme_min_buffered, geom_theme_od
        )
        geom_theme_od_min_clipped_plus_buffered = buffer_pos(
            geom_theme_od_clipped_min_buffered, self.relevant_distance
        )
        geom_theme_od_min_clipped_plus_buffered_clipped = safe_intersection(
            geom_theme_od_min_clipped_plus_buffered, geom_theme_od
        )
        geom_thematic_od = geom_theme_od_min_clipped_plus_buffered_clipped
        return geom_thematic_od

    def _od_snap_all_side(
        self,
        geometry,
        relevant_difference_array,
        relevant_intersection_array,
    ):
        geom_thematic_buffered = make_valid(
            buffer_pos(geometry, BUFFER_MULTIPLICATION_FACTOR * self.relevant_distance)
        )
        clip_ref_thematic_buffered = safe_intersection(
            self._get_reference_union(), geom_thematic_buffered
        )
        geom_reference = safe_difference(
            geom_thematic_buffered, clip_ref_thematic_buffered
        )  # beide OD-stukken worden SNAPPED toegevoegd
        geom_thematic = geometry
        if geom_reference.is_empty or geom_reference is None:
            geom_thematic_od = geom_reference
        else:
            geom_intersection = safe_intersection(geom_reference, geom_thematic)
            (
                geom_thematic_od,
                geom_relevant_intersection,
                geom_relevant_diff,
            ) = calculate_geom_by_intersection_and_reference(
                geom_intersection,
                geom_reference,
                True,
                self.relevant_distance / 2,
                self.threshold_overlap_percentage,
            )
            relevant_intersection_array = self.add_multi_polygons_from_geom_to_array(
                geom_relevant_intersection, relevant_intersection_array
            )
            relevant_difference_array = self.add_multi_polygons_from_geom_to_array(
                geom_relevant_diff, relevant_difference_array
            )
        return geom_thematic_od, relevant_difference_array, relevant_intersection_array

    # def _snap_geom_to_reference(self, geom_input, geom_reference, relevant_distance):
    # Deze functie werkt niet correct met Shapely. Deze vermijdt dat polygonen
    # collapsen als alles bijeen wordt gesnapt, wat we in sommige gevallen
    # effectief willen.
    # return snap(geom_input, geom_reference, relevant_distance)

    def _get_reference_union(self):
        if self.reference_union is None:
            self.reference_union = unary_union(list(self.dict_reference.values()))
        return self.reference_union

    def _postprocess_preresult(self, preresult, geom_thematic) -> ProcessResult:
        """
        Postprocess the preresult with the following actions to create the final result
        *Corrections for areas that differ more than the relevant distance
        *slivers
        *Inner holes (donuts) /multipolygons
        *validity
        *Circles (Polspy-popper)
        *Null/Empty-values

        Args:
            preresult (list): An existing list with all the elements of the preresult
            geom_thematic (BaseGeometry): The input geometry

        Returns:
            ProcessResult: A dictionary containing the resulting output geometries:

            *   result (BaseGeometry): The resulting output geometry
            *   result_diff (BaseGeometry): The resulting difference output geometry
            *   result_diff_plus (BaseGeometry): The resulting positive difference output
                geometry
            *   result_diff_min (BaseGeometry): The resulting negative difference output
                geometry
        """
        # Process array
        result = []
        geom_preresult = make_valid(unary_union(preresult))
        geom_thematic = make_valid(geom_thematic)

        # Corrections for areas that differ more than the relevant distance
        geom_thematic_dissolved = buffer_pos(
            buffer_neg(
                buffer_pos(geom_preresult, CORR_DISTANCE),
                2 * CORR_DISTANCE,
            ),
            CORR_DISTANCE,
        )
        # geom_symdiff = self._safe_symmetric_difference(geom_thematic,
        # geom_thematic_dissolved)
        geom_diff_add = safe_difference(geom_thematic, geom_thematic_dissolved)
        geom_diff_delete = safe_difference(geom_thematic_dissolved, geom_thematic)
        geom_diff_removed = safe_difference(
            geom_thematic_dissolved,
            safe_intersection(
                geom_diff_delete,
                buffer_neg_pos(geom_diff_delete, self.buffer_distance()),
            ),
        )
        geom_diff_removed_added = safe_union(
            geom_diff_removed,
            safe_intersection(
                geom_diff_add,
                buffer_neg_pos(geom_diff_add, self.buffer_distance()),
            ),
        )
        geom_thematic_preresult = buffer_pos(
            buffer_neg(
                buffer_pos(geom_diff_removed_added, CORR_DISTANCE),
                2 * CORR_DISTANCE,
            ),
            CORR_DISTANCE,
        )
        # Correction for Inner holes(donuts) / multipolygons
        # fill and remove gaps
        geom_thematic_cleaned_holes = geom_thematic_preresult
        ix_part = 1
        for part in get_parts(geom_thematic_preresult):
            exterior_ring = get_exterior_ring(part)
            exterior_polygon = polygons([exterior_ring])[0]
            empty_buffered_exterior_polygon = buffer_neg(
                exterior_polygon, self.buffer_distance()
            ).is_empty
            if (
                ix_part > 1
                and empty_buffered_exterior_polygon
                and not exterior_polygon.is_empty
            ):
                geom_thematic_cleaned_holes = safe_difference(
                    geom_thematic_cleaned_holes, exterior_polygon
                )
            num_interior_rings = get_num_interior_rings(part)
            if num_interior_rings > 0:
                ix = 0
                while ix < num_interior_rings:
                    interior_ring = get_interior_ring(part, ix)
                    interior_polygon = polygons([interior_ring])[0]

                    empty_buffered_interior_ring = buffer_neg(
                        interior_polygon, self.buffer_distance()
                    ).is_empty
                    if empty_buffered_interior_ring:
                        geom_thematic_cleaned_holes = safe_union(
                            geom_thematic_cleaned_holes, interior_polygon
                        )
                    ix = ix + 1
            ix_part = ix_part + 1

        geom_thematic_result = buffer_pos(
            buffer_neg(
                buffer_pos(geom_thematic_cleaned_holes, CORR_DISTANCE),
                2 * CORR_DISTANCE,
            ),
            CORR_DISTANCE,
        )
        geom_thematic_result = make_valid(remove_repeated_points(geom_thematic_result))
        # Correction for circles

        # calculate ratio to see if it is a circle, and keep the original geometry if a
        # circle: (Polsby-popper score)
        if not (geom_thematic.is_empty or geom_thematic is None):
            if (
                4 * pi * (geom_thematic.area / (geom_thematic.length**2))
                > THRESHOLD_CIRCLE_RATIO
            ):
                self.logger.feedback_warning(
                    "Circle: -->resulting geometry = original geometry"
                )
                geom_thematic_result = geom_thematic
        # Correction for empty preresults
        # TODO when a preresult gives an empty response, we better give an exception/warning instead of the original geometry
        if geom_thematic_result.is_empty or geom_thematic_result is None:
            self.logger.feedback_warning(
                "Empty result: -->resulting geometry = original geometry"
            )
            geom_thematic_result = geom_thematic

        # group all initial multipolygons into a new resulting dictionary
        result.append(geom_thematic_result)

        # create all resulting geometries
        geom_thematic_result = make_valid(unary_union(result))
        geom_result_diff = buffer_pos(
            buffer_neg(
                safe_symmetric_difference(geom_thematic_result, geom_thematic),
                CORR_DISTANCE,
            ),
            CORR_DISTANCE,
        )
        geom_result_diff_plus = safe_difference(geom_thematic_result, geom_thematic)
        geom_result_diff_min = safe_difference(geom_thematic, geom_thematic_result)

        return {
            "result": geom_thematic_result,
            "result_diff": geom_result_diff,
            "result_diff_plus": geom_result_diff_plus,
            "result_diff_min": geom_result_diff_min,
        }

    @staticmethod
    def add_multi_polygons_from_geom_to_array(geom: BaseGeometry, array):
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

    def load_reference_data(self, loader: Loader):
        self.dict_reference = loader.load_data()
        self._prepare_reference_data()

    def load_thematic_data(self, loader: Loader):
        self.dict_thematic = loader.load_data()

    # Deprecated loader methods
    def load_thematic_data_geojson(self, thematic_input, name_thematic_id):
        logging.warning("deprecated method, use load_thematic_data instead")
        loader = GeoJsonLoader(thematic_input, name_thematic_id)
        self.load_thematic_data(loader)

    def load_thematic_data_file(self, path_to_file, name_thematic_id):
        logging.warning("deprecated method, use load_thematic_data instead")
        loader = GeoJsonFileLoader(path_to_file, name_thematic_id)
        self.load_thematic_data(loader)

    def load_thematic_data_dict(self, dict_theme):
        logging.warning("deprecated method, use load_thematic_data instead")
        loader = DictLoader(dict_theme)
        self.load_thematic_data(loader)

    def load_thematic_data_url(self, url, name_thematic_id):
        logging.warning("deprecated method, use load_thematic_data instead")
        loader = GeoJsonUrlLoader(url, name_thematic_id)
        self.load_thematic_data(loader)

    def load_reference_data_dict(self, dict_ref):
        logging.warning("deprecated method, use load_reference_data instead")
        loader = DictLoader(dict_ref)
        self.load_reference_data(loader)

    def load_reference_data_geojson(self, reference_input, name_reference_id):
        logging.warning("deprecated method, use load_reference_data instead")
        loader = GeoJsonLoader(reference_input, name_reference_id)
        self.load_reference_data(loader)

    def load_reference_data_file(self, path_to_file, name_reference_id):
        logging.warning("deprecated method, use load_reference_data instead")
        loader = GeoJsonFileLoader(path_to_file, name_reference_id)
        self.load_reference_data(loader)

    def load_reference_data_url(self, url, name_reference_id):
        logging.warning("deprecated method, use load_reference_data instead")
        loader = GeoJsonUrlLoader(url, name_reference_id)
        self.load_reference_data(loader)

    def load_reference_data_grb_actual(self, *, grb_type=GRBType.ADP, partition=0):
        logging.warning("deprecated method, use load_reference_data instead")
        loader = GRBActualLoader(grb_type, partition, self)
        self.load_reference_data(loader)
