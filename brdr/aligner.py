import json
import logging
import os
from collections import defaultdict
from datetime import datetime
from math import pi
from typing import Iterable

import numpy as np
from shapely import GeometryCollection
from shapely import Polygon
from shapely import STRtree
from shapely import get_parts
from shapely import make_valid
from shapely import remove_repeated_points
from shapely import to_geojson
from shapely import unary_union
from shapely.geometry.base import BaseGeometry

from brdr import __version__
from brdr.constants import (
    BUFFER_MULTIPLICATION_FACTOR,
    LAST_VERSION_DATE,
    VERSION_DATE,
    DATE_FORMAT,
    THRESHOLD_EXCLUSION_PERCENTAGE,
    THRESHOLD_EXCLUSION_AREA,
    FORMULA_FIELD_NAME,
    EVALUATION_FIELD_NAME,
)
from brdr.constants import CORR_DISTANCE
from brdr.constants import DEFAULT_CRS
from brdr.constants import THRESHOLD_CIRCLE_RATIO
from brdr.enums import (
    OpenbaarDomeinStrategy,
    Evaluation,
    AlignerResultType,
    AlignerInputType,
)
from brdr.geometry_utils import buffer_neg
from brdr.geometry_utils import buffer_neg_pos
from brdr.geometry_utils import buffer_pos
from brdr.geometry_utils import fill_and_remove_gaps
from brdr.geometry_utils import safe_difference
from brdr.geometry_utils import safe_intersection
from brdr.geometry_utils import safe_symmetric_difference
from brdr.geometry_utils import safe_union
from brdr.loader import Loader
from brdr.logger import Logger
from brdr.typings import ProcessResult
from brdr.utils import diffs_from_dict_series, multipolygons_to_singles
from brdr.utils import geojson_from_dict
from brdr.utils import get_breakpoints_zerostreak
from brdr.utils import get_series_geojson_dict
from brdr.utils import merge_process_results
from brdr.utils import write_geojson


###################


class Aligner:
    """
    This class is used to compare the thematic data with the reference data.
    The reference data can be loaded in different ways, for example by using the GRB
    data.
    The thematic data can be loaded by using a geojson file.
    The class can be used to compare the thematic data with the reference data.
    """

    def __init__(
        self,
        *,
        feedback=None,
        relevant_distance=1,
        relevant_distances=np.arange(0, 200, 10, dtype=int) / 100,
        threshold_overlap_percentage=50,
        od_strategy=OpenbaarDomeinStrategy.SNAP_SINGLE_SIDE,
        crs=DEFAULT_CRS,
        area_limit=None,
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
            od_strategy (int, optional): Determines how the algorithm deals with parts
                of the geometry that are not on the
                reference (default 1: SNAP_SINGLE_SIDE)
            crs (str, optional): Coordinate Reference System (CRS) of the data.
                (default EPSG:31370)
            area_limit (int, optional): Maximum area for processing. (default 100000)


        """
        self.logger = Logger(feedback)
        self.relevant_distance = relevant_distance
        self.relevant_distances = relevant_distances
        self.od_strategy = od_strategy
        self.threshold_overlap_percentage = threshold_overlap_percentage
        self.area_limit = area_limit

        # PROCESSING DEFAULTS
        # thematic
        # name of the identifier-field of the thematic data (id has to be unique)
        self.name_thematic_id = "theme_identifier"
        # dictionary to store all thematic geometries to handle
        self.dict_thematic: dict[str, BaseGeometry] = {}
        # dictionary to store properties of the reference-features (optional)
        self.dict_thematic_properties: dict[str, dict] = {}
        # Dict to store source-information of the thematic dictionary
        self.dict_thematic_source: dict[str, str] = {}
        # dictionary to store all unioned thematic geometries
        self.thematic_union = None

        # reference

        # name of the identifier-field of the reference data (id has to be unique,f.e
        # CAPAKEY for GRB-parcels)
        self.name_reference_id = "ref_identifier"
        # dictionary to store all reference geometries
        self.dict_reference: dict[str, BaseGeometry] = {}
        # dictionary to store properties of the reference-features (optional)
        self.dict_reference_properties: dict[str, dict] = {}
        # Dict to store source-information of the reference dictionary
        self.dict_reference_source: dict[str, str] = {}
        # to save a unioned geometry of all reference polygons; needed for calculation
        # in most OD-strategies
        self.reference_union = None

        # results

        # output-dictionaries (all results of process()), grouped by theme_id and relevant_distance
        self.dict_processresults: dict[str, dict[float, ProcessResult]] = {}
        # dictionary with the 'predicted' results, grouped by theme_id and relevant_distance
        self.dict_predictions: dict[str, dict[float, ProcessResult]] = {}

        # Coordinate reference system
        # thematic geometries and reference geometries are assumed to be in the same CRS
        # before loading into the Aligner. No CRS-transformation will be performed.
        # When loading data, CRS is expected to be a projected CRS with units in 'meter
        # (m)'.
        # Default EPSG:31370 (Lambert72), alternative: EPSG:3812 (Lambert2008)
        self.CRS = crs
        # this parameter is used to treat multipolygon as single polygons. So polygons
        # with ID splitter are separately evaluated and merged on result.
        self.multi_as_single_modus = True
        self.logger.feedback_info("Aligner initialized")

    def buffer_distance(self):
        return self.relevant_distance / 2

    ##########LOADERS##########################
    ###########################################

    def load_thematic_data(self, loader: Loader):
        self.dict_thematic, self.dict_thematic_properties, self.dict_thematic_source = (
            loader.load_data()
        )

    def load_reference_data(self, loader: Loader):
        (
            self.dict_reference,
            self.dict_reference_properties,
            self.dict_reference_source,
        ) = loader.load_data()
        self._prepare_reference_data()

    ##########PROCESSORS#######################
    ###########################################

    def process_geometry(
        self,
        input_geometry: BaseGeometry,
        relevant_distance: float = 1,
        od_strategy=OpenbaarDomeinStrategy.SNAP_SINGLE_SIDE,
        threshold_overlap_percentage=50,
    ) -> ProcessResult:
        """
        method to align a geometry to the reference layer

        Args:
            input_geometry (BaseGeometry): The input geometric object.
            relevant_distance
            od_strategy
            threshold_overlap_percentage (int): The buffer distance (positive
                or negative).

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
        Notes:
            -
        Example:
        """
        if self.area_limit and input_geometry.area > self.area_limit:
            message = "The input geometry is too large to process."
            raise ValueError(message)

        self.logger.feedback_debug("process geometry")
        self.relevant_distance = relevant_distance
        self.od_strategy = od_strategy
        self.threshold_overlap_percentage = threshold_overlap_percentage

        # combine all parts of the input geometry to one polygon
        input_geometry = unary_union(get_parts(input_geometry))

        # array with all relevant parts of a thematic geometry; initial empty Polygon
        (
            geometry,
            preresult,
            relevant_intersection_array,
            relevant_diff_array,
        ) = self._calculate_intersection_between_geometry_and_od(input_geometry)
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
            ) = _calculate_geom_by_intersection_and_reference(
                geom_intersection,
                geom_reference,
                False,
                self.relevant_distance / 2,
                self.threshold_overlap_percentage,
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

        # make a unary union for each key value in the result dict
        for key in ProcessResult.__annotations__:
            geometry = result_dict.get(key, Polygon())  # noqa
            if not geometry.is_empty:
                geometry = unary_union(geometry)
            result_dict[key] = geometry  # noqa

        return result_dict

    def process(
        self,
        relevant_distances: Iterable[float] = None,
        relevant_distance=1,
        od_strategy=OpenbaarDomeinStrategy.SNAP_SINGLE_SIDE,
        threshold_overlap_percentage=50,
    ) -> dict[str, dict[float, ProcessResult]]:
        """
        Calculates the resulting dictionaries for thematic data based on a series of
            relevant distances.

        Args:
            relevant_distances (Iterable[float]): A series of relevant distances
                (in meters) to process
            od_strategy (int, optional): The strategy for overlap detection.
                Defaults to 1.
            threshold_overlap_percentage (float, optional): The threshold percentage for
                considering full overlap. Defaults to 50.

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
            self.relevant_distance = relevant_distance
        self.relevant_distances = relevant_distances
        self.od_strategy = od_strategy
        self.threshold_overlap_percentage = threshold_overlap_percentage
        self.logger.feedback_debug("Process series" + str(self.relevant_distances))
        dict_series = {}
        dict_thematic = self.dict_thematic

        if self.multi_as_single_modus:
            dict_thematic = multipolygons_to_singles(dict_thematic)

        for key, geometry in dict_thematic.items():
            self.logger.feedback_info(
                f"thematic id {str(key)} processed with relevant distances (m) [{str(self.relevant_distances)}]"
            )
            dict_series[key] = {}
            for relevant_distance in self.relevant_distances:
                try:
                    self.relevant_distance = relevant_distance
                    processed_result = self.process_geometry(
                        geometry,
                        self.relevant_distance,
                        od_strategy,
                        threshold_overlap_percentage,
                    )
                except ValueError as e:
                    self.logger.feedback_warning(str(e))

                dict_series[key][self.relevant_distance] = processed_result

        if self.multi_as_single_modus:
            dict_series = merge_process_results(dict_series)

        self.logger.feedback_info(
            "End of processing series: " + str(self.relevant_distances)
        )
        self.dict_processresults = dict_series

        return self.dict_processresults

    # def process(
    #     self,
    #     relevant_distance=1,
    #     od_strategy=OpenbaarDomeinStrategy.SNAP_SINGLE_SIDE,
    #     threshold_overlap_percentage=50,
    # ) -> dict[str, dict[float, ProcessResult]]:
    #     """
    #     Aligns a thematic dictionary of geometries to the reference layer based on
    #     specified parameters. - method to align a thematic dictionary to the reference
    #     layer
    #
    #     Args:
    #         relevant_distance (float, optional): The relevant distance (in meters) for
    #             processing. Defaults to 1.
    #         od_strategy (int, optional): The strategy for overlap detection.
    #             Defaults to 1.
    #         threshold_overlap_percentage (float, optional): The threshold percentage for
    #             considering full overlap. Defaults to 50.
    #
    #     Returns:
    #         dict: A dict containing processed data for each thematic key:
    #             - result: Aligned thematic data.
    #             - result_diff: global differences between thematic data and reference
    #               data.
    #             - result_diff_plus: Positive differences.
    #             - result_diff_min: Negative differences.
    #             - relevant_intersection: relevant intersections.
    #             - relevant_diff: relevant differences.
    #
    #     """
    #     self.relevant_distance=relevant_distance
    #     self.dict_result = self.process(relevant_distances=[self.relevant_distance],
    #                                     od_strategy=od_strategy,
    #                                     threshold_overlap_percentage=threshold_overlap_percentage)
    #     return self.dict_result

    def predictor(
        self,
        relevant_distances=np.arange(0, 300, 10, dtype=int) / 100,
        od_strategy=OpenbaarDomeinStrategy.SNAP_SINGLE_SIDE,
        threshold_overlap_percentage=50,
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

        5. **Filter Results:**
            - The function might further filter the predicted results for each thematic
              element based on the element key (using `filter_resulting_series_by_key`).

        Args:
            relevant_distances (np.ndarray, optional): A NumPy array of distances to
              be analyzed. Defaults to np.arange(0.1, 5.05, 0.1).
            od_strategy (OpenbaarDomeinStrategy, optional): A strategy for handling
              open data in the processing (implementation specific). Defaults to
             OpenbaarDomeinStrategy.SNAP_SINGLE_SIDE.
            threshold_overlap_percentage (int, optional): A percentage threshold for
              considering full overlap in the processing (implementation specific).
             Defaults to 50.

        Returns:
            dict: A dictionary containing predicted interesting distances for each
            thematic element.
                - Keys: Thematic element identifiers from `self.dict_thematic`.
                - Values: Dictionaries with the following structure for each
                   thematic element:
                    - Keys: Distances identified as interesting (breakpoints or
                    zero-streaks).
                    - Values: dicts containing results (likely specific to
                    your implementation) from the distance series for the
                    corresponding distance.

        Logs:
            - Debug logs the thematic element key being processed.
        """
        dict_predictions = defaultdict(dict)
        dict_series = self.process(
            relevant_distances=relevant_distances,
            od_strategy=od_strategy,
            threshold_overlap_percentage=threshold_overlap_percentage,
        )
        dict_thematic = self.dict_thematic

        diffs_dict = diffs_from_dict_series(dict_series, dict_thematic)

        for theme_id, diffs in diffs_dict.items():
            if len(diffs) != len(relevant_distances):
                logging.warning(
                    f"Number of computed diffs for thematic element {theme_id} does "
                    f"not match the number of relevant distances."
                )
                continue
            diff_values = list(diffs.values())
            breakpoints, zero_streaks = get_breakpoints_zerostreak(
                relevant_distances, diff_values
            )
            logging.debug(str(theme_id))
            if len(zero_streaks) == 0:
                dict_predictions[theme_id][relevant_distances[0]] = dict_series[
                    theme_id
                ][relevant_distances[0]]
                logging.info("No zero-streaks found for: " + str(theme_id))
            for zs in zero_streaks:
                dict_predictions[theme_id][zs[0]] = dict_series[theme_id][zs[0]]

        # Check if the predicted reldists are unique (and remove duplicated predictions
        dict_predictions_unique = defaultdict(dict)
        for theme_id, dist_results in dict_predictions.items():
            dict_predictions_unique[theme_id] = {}
            predicted_geoms_for_theme_id = []
            for rel_dist, processresults in dist_results.items():
                predicted_geom = processresults["result"]
                if not _equal_geom_in_array(
                    predicted_geom, predicted_geoms_for_theme_id
                ):
                    dict_predictions_unique[theme_id][rel_dist] = processresults
                    predicted_geoms_for_theme_id.append(processresults["result"])
                else:
                    self.logger.feedback_info(
                        f"Duplicate prediction found for key {theme_id} at distance {rel_dist}: Prediction excluded"
                    )

        self.dict_predictions = dict_predictions_unique

        return (
            dict_series,
            self.dict_predictions,
            diffs_dict,
        )

    def compare(
        self,
        threshold_area=5,
        threshold_percentage=1,
        dict_unchanged=None,
    ):
        """
        Compares input-geometries (with formula) and evaluates these geometries: An attribute is added to evaluate and decide if new
        proposals can be used
        """
        dict_series, dict_predictions, diffs = self.predictor(self.relevant_distances)
        if dict_unchanged is None:
            dict_unchanged = {}
        theme_ids = list(dict_series.keys())
        dict_evaluated_result = {}
        prop_dictionary = {}
        # Fill the dictionary-structure with empty values
        for theme_id in theme_ids:
            dict_evaluated_result[theme_id] = {}
            prop_dictionary[theme_id] = {}
            for dist in dict_series[theme_id].keys():
                prop_dictionary[theme_id][dist] = {}
        for theme_id in dict_unchanged.keys():
            prop_dictionary[theme_id] = {}

        for theme_id, dict_results in dict_predictions.items():
            equality = False
            for dist in sorted(dict_results.keys()):
                if equality:
                    break
                geomresult = dict_results[dist]["result"]
                actual_formula = self.get_brdr_formula(geomresult)
                prop_dictionary[theme_id][dist][FORMULA_FIELD_NAME] = json.dumps(
                    actual_formula
                )
                base_formula = None
                if (
                    theme_id in self.dict_thematic_properties
                    and FORMULA_FIELD_NAME in self.dict_thematic_properties[theme_id]
                ):
                    base_formula = self.dict_thematic_properties[theme_id][
                        FORMULA_FIELD_NAME
                    ]
                equality, prop = _check_equality(
                    base_formula,
                    actual_formula,
                    threshold_area,
                    threshold_percentage,
                )
                if equality:
                    dict_evaluated_result[theme_id][dist] = dict_predictions[theme_id][
                        dist
                    ]
                    prop_dictionary[theme_id][dist][EVALUATION_FIELD_NAME] = prop
                    break

        evaluated_theme_ids = [
            theme_id for theme_id, value in dict_evaluated_result.items() if value != {}
        ]

        # fill where no equality is found/ The biggest predicted distance is returned as
        # proposal
        for theme_id in theme_ids:
            if theme_id not in evaluated_theme_ids:
                if len(dict_predictions[theme_id].keys()) == 0:
                    result = dict_series[theme_id][0]
                    dict_evaluated_result[theme_id][0] = result
                    prop_dictionary[theme_id][0][FORMULA_FIELD_NAME] = json.dumps(
                        self.get_brdr_formula(result["result"])
                    )
                    prop_dictionary[theme_id][0][
                        EVALUATION_FIELD_NAME
                    ] = Evaluation.NO_PREDICTION_5
                    continue
                # Add all predicted features so they can be manually checked
                for dist in dict_predictions[theme_id].keys():
                    predicted_resultset = dict_predictions[theme_id][dist]
                    dict_evaluated_result[theme_id][dist] = predicted_resultset
                    prop_dictionary[theme_id][dist][FORMULA_FIELD_NAME] = json.dumps(
                        self.get_brdr_formula(predicted_resultset["result"])
                    )
                    prop_dictionary[theme_id][dist][
                        EVALUATION_FIELD_NAME
                    ] = Evaluation.TO_CHECK_4

        for theme_id, geom in dict_unchanged.items():
            prop_dictionary[theme_id] = {
                0: {
                    "result": geom,
                    EVALUATION_FIELD_NAME: Evaluation.NO_CHANGE_6,
                    FORMULA_FIELD_NAME: json.dumps(self.get_brdr_formula(geom)),
                }
            }
        return dict_evaluated_result, prop_dictionary

    def get_brdr_formula(self, geometry: BaseGeometry, with_geom=False):
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
        dict_formula = {
            "alignment_date": datetime.now().strftime(DATE_FORMAT),
            "brdr_version": str(__version__),
            "reference_source": self.dict_reference_source,
            "full": True,
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

            perc = round(geom_intersection.area * 100 / geom_reference.area, 2)
            if perc < 0.01:
                continue
            # Add a last_version_date if available in properties
            if (
                key_ref in self.dict_reference_properties
                and VERSION_DATE in self.dict_reference_properties[key_ref]
            ):
                version_date = self.dict_reference_properties[key_ref][VERSION_DATE]
                if last_version_date is None and version_date is not None:
                    last_version_date = version_date
                if version_date is not None and version_date > last_version_date:
                    last_version_date = version_date

            if perc > 99.99:
                full = True
                area = round(geom_reference.area, 2)
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
                safe_difference(geometry, make_valid(unary_union(intersected))),
                CORR_DISTANCE,
            ),
            CORR_DISTANCE,
        )
        if geom_od is not None:
            area_od = round(geom_od.area, 2)
            if area_od > 0:
                dict_formula["reference_od"] = {"area": area_od}
                if with_geom:
                    dict_formula["reference_od"]["geometry"] = to_geojson(geom_od)
        self.logger.feedback_debug(str(dict_formula))
        return dict_formula

    ##########EXPORTERS########################
    ###########################################

    def get_results_as_geojson(
        self, resulttype=AlignerResultType.PROCESSRESULTS, formula=False
    ):
        """
        get a geojson of  a dictionary containing the resulting geometries for all
            'serial' relevant distances. If no dict_series is given, the dict_result returned.
        Optional: The descriptive formula is added as an attribute to the result"""

        if resulttype == AlignerResultType.PROCESSRESULTS:
            dict_series = self.dict_processresults
        elif resulttype == AlignerResultType.PREDICTIONS:
            dict_series = self.dict_predictions
        else:
            raise (ValueError, "AlignerResultType unknown")
        if dict_series is None or dict_series == {}:
            self.logger.feedback_warning(
                "Empty results: No calculated results to export."
            )
            return {}

        prop_dictionary = defaultdict(dict)

        for theme_id, results_dict in dict_series.items():
            for relevant_distance, process_results in results_dict.items():
                if formula:
                    result = process_results["result"]
                    formula = self.get_brdr_formula(result)
                    prop_dictionary[theme_id][relevant_distance] = {
                        FORMULA_FIELD_NAME: json.dumps(formula)
                    }

        return get_series_geojson_dict(
            dict_series,
            crs=self.CRS,
            id_field=self.name_thematic_id,
            series_prop_dict=prop_dictionary,
        )

    def get_input_as_geojson(self, inputtype=AlignerInputType.REFERENCE):
        """
        get a geojson of the reference polygons
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
        dict_properties
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
        Exports analysis results as GeoJSON files.

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
        if self.thematic_union is None:
            self.thematic_union = make_valid(
                unary_union(list(self.dict_thematic.values()))
            )
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
        self.reference_items = np.array(list(self.dict_reference.keys()))
        # clear the reference_union, so it will be recalculated on request when needed
        self.reference_union = None
        return

    def _calculate_intersection_between_geometry_and_od(self, geometry):
        # Calculate the intersection between thematic and Openbaar Domein
        relevant_intersection_array = []
        relevant_difference_array = []
        geom_thematic_od = Polygon()

        if self.od_strategy == OpenbaarDomeinStrategy.EXCLUDE:
            # Completely exclude everything that is not on the reference layer
            self.logger.feedback_debug("OD-strategy EXCLUDE")
            # Remove from the thematic layer all parts that are not on the reference
            # layer
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
            ) = self._od_snap_all_side(geometry)
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
            ) = self._od_snap_all_side(geometry)
            # This part is a copy of SNAP_FULL_AREA_SINGLE_SIDE
            geom_theme_od_min_clipped_plus_buffered_clipped = self._od_full_area(
                geometry
            )
            # UNION the calculation of  OD-SNAP_ALL_SIDE with FULL AREA of
            # OD-SNAP_FULL_AREA_SINGLE_SIDE
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
            geom_thematic_od = Polygon()
            pass

        # ADD THEMATIC_OD
        preresult = self._add_multi_polygons_from_geom_to_array(geom_thematic_od, [])
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

    def _od_snap_all_side(self, geometry):
        relevant_difference_array = []
        relevant_intersection_array = []
        geom_thematic_buffered = make_valid(
            buffer_pos(geometry, BUFFER_MULTIPLICATION_FACTOR * self.relevant_distance)
        )
        clip_ref_thematic_buffered = safe_intersection(
            self._get_reference_union(), geom_thematic_buffered
        )
        geom_reference = safe_difference(
            geom_thematic_buffered, clip_ref_thematic_buffered
        )  # Both OD-parts are SNAPPED added
        geom_thematic = geometry
        if geom_reference.is_empty or geom_reference is None:
            geom_thematic_od = geom_reference
        else:
            geom_intersection = safe_intersection(geom_reference, geom_thematic)
            (
                geom_thematic_od,
                geom_relevant_intersection,
                geom_relevant_diff,
            ) = _calculate_geom_by_intersection_and_reference(
                geom_intersection,
                geom_reference,
                True,
                self.relevant_distance / 2,
                self.threshold_overlap_percentage,
            )
            relevant_intersection_array = self._add_multi_polygons_from_geom_to_array(
                geom_relevant_intersection, []
            )
            relevant_difference_array = self._add_multi_polygons_from_geom_to_array(
                geom_relevant_diff, []
            )
        return geom_thematic_od, relevant_difference_array, relevant_intersection_array

    def _get_reference_union(self):
        if self.reference_union is None:
            self.reference_union = make_valid(
                unary_union(list(self.dict_reference.values()))
            )
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
            *   result_diff_plus (BaseGeometry): The resulting positive difference
                output geometry
            *   result_diff_min (BaseGeometry): The resulting negative difference output
                geometry
        """
        # Process array
        result = []
        geom_preresult = make_valid(unary_union(preresult))
        geom_thematic = make_valid(geom_thematic)

        if not (geom_thematic is None or geom_thematic.is_empty):
            # Correction for circles
            # calculate ratio to see if it is a circle, and keep the original geometry
            #  if a circle: (Polsby-popper score)
            if (
                4 * pi * (geom_thematic.area / (geom_thematic.length**2))
                > THRESHOLD_CIRCLE_RATIO
            ):
                self.logger.feedback_debug(
                    "Circle: -->resulting geometry = original geometry"
                )
                return {"result": geom_thematic}

            # Correction for unchanged geometries
            if geom_preresult == geom_thematic:
                return {"result": geom_thematic}

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
        geom_thematic_cleaned_holes = fill_and_remove_gaps(
            geom_thematic_preresult, self.buffer_distance()
        )
        geom_thematic_result = buffer_pos(
            buffer_neg(
                buffer_pos(geom_thematic_cleaned_holes, CORR_DISTANCE),
                2 * CORR_DISTANCE,
            ),
            CORR_DISTANCE,
        )
        geom_thematic_result = make_valid(remove_repeated_points(geom_thematic_result))

        # Correction for empty preresults
        if geom_thematic_result.is_empty or geom_thematic_result is None:
            self.logger.feedback_warning(
                "Empty result: -->resulting geometry = empty geometry"
            )
            # geom_thematic_result = geom_thematic
            geom_thematic_result = Polygon()

        # group all initial multipolygons into a new resulting dictionary
        result.append(geom_thematic_result)

        # create all resulting geometries
        geom_thematic_result = make_valid(unary_union(result))

        # negative and positive buffer is added to the difference-calculations, to
        # remove 'very small' differences (smaller than the correction distance)
        geom_result_diff = buffer_pos(
            buffer_neg(
                safe_symmetric_difference(geom_thematic_result, geom_thematic),
                CORR_DISTANCE,
            ),
            CORR_DISTANCE,
        )
        geom_result_diff_plus = buffer_pos(
            buffer_neg(
                safe_difference(geom_thematic_result, geom_thematic),
                CORR_DISTANCE,
            ),
            CORR_DISTANCE,
        )
        geom_result_diff_min = buffer_pos(
            buffer_neg(
                safe_difference(geom_thematic, geom_thematic_result),
                CORR_DISTANCE,
            ),
            CORR_DISTANCE,
        )
        # geom_result_diff_plus = safe_difference(geom_thematic_result, geom_thematic)
        # geom_result_diff_min = safe_difference(geom_thematic, geom_thematic_result)

        return {
            "result": geom_thematic_result,
            "result_diff": geom_result_diff,
            "result_diff_plus": geom_result_diff_plus,
            "result_diff_min": geom_result_diff_min,
        }

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
    is_openbaar_domein,
    buffer_distance,
    threshold_overlap_percentage,
    threshold_exclusion_percentage=THRESHOLD_EXCLUSION_PERCENTAGE,
    threshold_exclusion_area=THRESHOLD_EXCLUSION_AREA,
):
    """
    Calculates the geometry based on intersection and reference geometries.

    Args:
        geom_intersection (BaseGeometry): The intersection geometry.
        geom_reference (BaseGeometry): The reference geometry.
        is_openbaar_domein (bool): A flag indicating whether it's a public domain
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

    if geom_reference.area == 0:
        overlap = 100

    else:
        overlap = geom_intersection.area * 100 / geom_reference.area

    if (
        overlap < threshold_exclusion_percentage
        or geom_intersection.area < threshold_exclusion_area
    ):
        return Polygon(), Polygon(), Polygon()

    geom_difference = safe_difference(geom_reference, geom_intersection)
    geom_relevant_intersection = buffer_neg(geom_intersection, buffer_distance)
    geom_relevant_difference = buffer_neg(geom_difference, buffer_distance)
    if (
        not geom_relevant_intersection.is_empty
        and not geom_relevant_difference.is_empty
    ):
        # relevant intersection and relevant difference
        geom_x = safe_intersection(
            geom_reference,
            safe_difference(
                geom_reference,
                safe_intersection(
                    geom_difference,
                    buffer_neg_pos(geom_difference, buffer_distance),
                ),
            ),
        )
        geom = safe_intersection(
            geom_x,
            buffer_pos(
                buffer_neg_pos(geom_x, buffer_distance),
                buffer_distance,
            ),
        )
        # when calculating for OD, we create a 'virtual parcel'. When calculating this
        # virtual parcel, it is buffered to take outer boundaries into account.
        # This results in a side effect that there are extra non-logical parts included
        # in the result. The function below tries to exclude these non-logical parts.
        # see eo_id 206363 with relevant distance=0.2m and SNAP_ALL_SIDE
        if is_openbaar_domein:
            geom = _get_relevant_polygons_from_geom(geom, buffer_distance)
    elif not geom_relevant_intersection.is_empty and geom_relevant_difference.is_empty:
        geom = geom_reference
    elif geom_relevant_intersection.is_empty and not geom_relevant_difference.is_empty:
        geom = geom_relevant_intersection  # (=empty geometry)
    else:
        if is_openbaar_domein:
            geom = geom_relevant_intersection  # (=empty geometry)
        # geom = snap_geom_to_reference (geom_intersection, geom_reference,
        # relevant_distance)
        elif threshold_overlap_percentage < 0:
            # if we take a value of -1, the original border will be used
            geom = geom_intersection
        elif overlap > threshold_overlap_percentage:
            geom = geom_reference
        else:
            geom = geom_relevant_intersection  # (=empty geometry)
    return geom, geom_relevant_intersection, geom_relevant_difference


def _get_relevant_polygons_from_geom(geometry: BaseGeometry, buffer_distance: float):
    """
    Get only the relevant parts (polygon) from a geometry.
    Points, Lines and Polygons smaller than relevant distance are excluded from the
    result
    """
    if not geometry or geometry.is_empty:
        # If the input geometry is empty or None, do nothing.
        return geometry
    else:
        geometry = make_valid(unary_union(geometry))
        # Create a GeometryCollection from the input geometry.
        geometry_collection = GeometryCollection(geometry)
        array = []
        for g in geometry_collection.geoms:
            # Ensure each sub-geometry is valid.
            g = make_valid(g)
            if str(g.geom_type) in ["Polygon", "MultiPolygon"]:
                relevant_geom = buffer_neg(g, buffer_distance)
                if relevant_geom is not None and not relevant_geom.is_empty:
                    array.append(g)
    return make_valid(unary_union(array))


def _equal_geom_in_array(geom, geom_array):
    """
    Check if a predicted geometry is equal to other predicted geometries in a list.
    Equality is defined as there is the symmetrical difference is smaller than the CORRECTION DISTANCE
    Returns True if one of the elements is equal, otherwise False
    """
    for g in geom_array:
        # if safe_equals(geom,g):
        if buffer_neg(safe_symmetric_difference(geom, g), CORR_DISTANCE).is_empty:
            return True
    return False


def _check_equality(
    base_formula, actual_formula, threshold_area=5, threshold_percentage=1
):
    """
    function that checks if 2 formulas are equal (True,False) and adds an Evaluation
    """
    if base_formula is None or actual_formula is None:
        return False, Evaluation.NO_PREDICTION_5
    od_alike = False
    if base_formula["reference_od"] is None and actual_formula["reference_od"] is None:
        od_alike = True
    elif base_formula["reference_od"] is None or actual_formula["reference_od"] is None:
        od_alike = False
    elif (
        abs(
            base_formula["reference_od"]["area"]
            - actual_formula["reference_od"]["area"]
        )
        * 100
        / base_formula["reference_od"]["area"]
    ) < threshold_percentage:
        od_alike = True

    if (
        base_formula["reference_features"].keys()
        == actual_formula["reference_features"].keys()
        and od_alike
    ):
        if base_formula["full"] and base_formula["full"]:
            return True, Evaluation.EQUALITY_FORMULA_GEOM_1

        equal_reference_features = True
        for key in base_formula["reference_features"].keys():
            if (
                (
                    base_formula["reference_features"][key]["full"]
                    == actual_formula["reference_features"][key]["full"]
                )
                or (
                    abs(
                        base_formula["reference_features"][key]["area"]
                        - actual_formula["reference_features"][key]["area"]
                    )
                    > threshold_area
                )
                or (
                    (
                        abs(
                            base_formula["reference_features"][key]["area"]
                            - actual_formula["reference_features"][key]["area"]
                        )
                        * 100
                        / base_formula["reference_features"][key]["area"]
                    )
                    > threshold_percentage
                )
            ):
                equal_reference_features = False
        if equal_reference_features:
            return True, Evaluation.EQUALITY_FORMULA_2
    if base_formula["full"] and base_formula["full"] and od_alike:
        return True, Evaluation.EQUALITY_GEOM_3
    return False, Evaluation.NO_PREDICTION_5
