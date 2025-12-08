import json
import os
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import wait
from datetime import datetime
from typing import Iterable

import numpy as np
from shapely import GeometryCollection
from shapely import STRtree
from shapely import make_valid
from shapely import to_geojson
from shapely.geometry.base import BaseGeometry

from brdr import __version__
from brdr.configs import ProcessorConfig
from brdr.constants import DATE_FORMAT
from brdr.constants import DEFAULT_CRS
from brdr.constants import DIFF_AREA_FIELD_NAME
from brdr.constants import DIFF_INDEX
from brdr.constants import DIFF_METRIC
from brdr.constants import DIFF_PERCENTAGE_FIELD_NAME
from brdr.constants import DIFF_PERC_INDEX
from brdr.constants import EQUAL_REFERENCE_FEATURES_FIELD_NAME
from brdr.constants import EVALUATION_FIELD_NAME
from brdr.constants import FORMULA_FIELD_NAME
from brdr.constants import FULL_ACTUAL_FIELD_NAME
from brdr.constants import FULL_BASE_FIELD_NAME
from brdr.constants import ID_REFERENCE_FIELD_NAME
from brdr.constants import ID_THEME_FIELD_NAME
from brdr.constants import LAST_VERSION_DATE
from brdr.constants import NR_CALCULATION_FIELD_NAME
from brdr.constants import OD_ALIKE_FIELD_NAME
from brdr.constants import PREDICTION_COUNT
from brdr.constants import PREDICTION_SCORE
from brdr.constants import RELEVANT_DISTANCE_DECIMALS
from brdr.constants import RELEVANT_DISTANCE_FIELD_NAME
from brdr.constants import REMARK_FIELD_NAME
from brdr.constants import STABILITY
from brdr.constants import VERSION_DATE
from brdr.constants import ZERO_STREAK
from brdr.enums import AlignerInputType
from brdr.enums import AlignerResultType
from brdr.enums import DiffMetric
from brdr.enums import Evaluation
from brdr.enums import FullStrategy
from brdr.geometry_utils import buffer_neg
from brdr.geometry_utils import buffer_pos
from brdr.geometry_utils import extract_points_lines_from_geometry
from brdr.geometry_utils import safe_difference
from brdr.geometry_utils import safe_intersection
from brdr.geometry_utils import safe_unary_union
from brdr.geometry_utils import to_crs
from brdr.loader import Loader
from brdr.logger import Logger
from brdr.processor import AlignerGeometryProcessor
from brdr.processor import BaseProcessor
from brdr.topo import dissolve_topo
from brdr.topo import generate_topo
from brdr.typings import ProcessResult
from brdr.typings import ThematicId
from brdr.utils import coverage_ratio
from brdr.utils import determine_stability
from brdr.utils import diff_from_processresult
from brdr.utils import diffs_from_dict_processresult
from brdr.utils import equal_geom_in_array
from brdr.utils import geojson_from_dict
from brdr.utils import get_dict_geojsons_from_series_dict
from brdr.utils import is_brdr_formula
from brdr.utils import merge_process_results
from brdr.utils import multi_to_singles
from brdr.utils import write_geojson


###################
class AlignerResult:
    def __init__(
        self,
        process_results: dict[ThematicId, dict[float, ProcessResult | None]],
    ):
        self.results = process_results

    def get_results(
        self,
        result_type:AlignerResultType=AlignerResultType.PROCESSRESULTS,
    ) -> dict[ThematicId, dict[float, ProcessResult | None]]:
        if result_type == AlignerResultType.PROCESSRESULTS:
            return self.results
        elif result_type == AlignerResultType.PREDICTIONS:
            # return all theme ids with only the relevant distances
            # where the ProcessResult has property PREDICTION_SCORE:
            return {
                theme_id: {
                    rd: process_result
                    for rd, process_result in results_dict.items()
                    if process_result is not None
                    and PREDICTION_SCORE in process_result["properties"]
                }
                for theme_id, results_dict in self.results.items()
            }
        elif result_type == AlignerResultType.EVALUATED_PREDICTIONS:
            # return all theme ids with only the relevant distances
            # where the ProcessResult has property EVALUATION_FIELD_NAME:
            return {
                theme_id: {
                    rd: process_result
                    for rd, process_result in results_dict.items()
                    if process_result is not None
                    and EVALUATION_FIELD_NAME in process_result["properties"]
                }
                for theme_id, results_dict in self.results.items()
            }
        else:
            raise ValueError("Unknown result type: " + str(result_type))

    def get_results_as_geojson(
        self,
        aligner,
        formula=False,
        attributes=False,
    ):
        """
        get a geojson of a dictionary containing the resulting geometries for all
            'serial' relevant distances. The resulttype can be chosen.
        formula (boolean, Optional): The descriptive formula is added as an attribute to the result
        attributes (boolean, Optional): The original attributes/properties are added to the result
        """

        if self.results is None or self.results == {}:
            raise ValueError("Empty results: No calculated results to export.")

        prop_dictionary = defaultdict(dict)

        for theme_id, results_dict in self.results.items():
            nr_calculations = len(results_dict)
            for relevant_distance, process_result in results_dict.items():
                properties = process_result["properties"]

                # Adding extra properties
                properties[ID_THEME_FIELD_NAME] = theme_id
                properties[NR_CALCULATION_FIELD_NAME] = nr_calculations
                properties[RELEVANT_DISTANCE_FIELD_NAME] = relevant_distance
                properties[DIFF_INDEX] = diff_from_processresult(
                    process_result,
                    aligner.dict_thematic[theme_id],
                    None,
                    DIFF_METRIC.CHANGES_AREA,
                )
                properties[DIFF_PERC_INDEX] = diff_from_processresult(
                    process_result,
                    aligner.dict_thematic[theme_id],
                    None,
                    DIFF_METRIC.CHANGES_PERCENTAGE,
                )
                prop_dictionary[theme_id][relevant_distance] = properties

                # Adding original attributes
                if attributes and theme_id in aligner.dict_thematic_properties.keys():
                    for attr, value in aligner.dict_thematic_properties[
                        theme_id
                    ].items():
                        prop_dictionary[theme_id][relevant_distance][attr] = value

                # Adding formula
                if (
                    formula
                ):  # and not (theme_id in prop_dictionary and relevant_distance in prop_dictionary[theme_id] and NEW_FORMULA_FIELD_NAME in prop_dictionary[theme_id][relevant_distance]):
                    result = process_result["result"]
                    formula = aligner.get_brdr_formula(result)
                    prop_dictionary[theme_id][relevant_distance][FORMULA_FIELD_NAME] = (
                        json.dumps(formula)
                    )
        return get_dict_geojsons_from_series_dict(
            self.results,
            crs=aligner.CRS,
            id_field=aligner.name_thematic_id,
            series_prop_dict=prop_dictionary,
        )

    def save_results(self, aligner, path, formula=True, result_type = AlignerResultType.PROCESSRESULTS):
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
            aligner=aligner,
            formula=formula,
        )
        for name, fc in fcs.items():
            write_geojson(
                os.path.join(path, result_type.value + "_" + name + ".geojson"), fc
            )


# TODO what about the Aligner-parameters; AlignerConfig-class?
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
        processor: BaseProcessor = None,
        crs=DEFAULT_CRS,
        multi_as_single_modus=True,
        preserve_topology=False,
        correction_distance=0.01,
        diff_metric=DIFF_METRIC,
        mitre_limit=10,
        max_workers=None,
    ):
        """
        Initializes the Aligner object

        Args:
            feedback (object, optional): Feedback object that can be added to show
                feedback in QGIS. Defaults to None.
            crs (str, optional): Coordinate Reference System (CRS) of the data.
                (default EPSG:31370)
            multi_as_single_modus (boolean, optional): Modus to handle multipolygons (Default=True):
                True: input-multipolygons will be split-up into single polygons and handled by the algorithm. After executing the algorithm, the results are merged together.
                False: Multipolygons are directly processed by the algorithm
            correction_distance (float, optional): Distance used in a pos_neg_buffer to remove slivers (technical correction) (Default= 0.01 = 1cm )
            mitre_limit (int, optional):buffer-parameter - The mitre ratio is the ratio of the distance from the corner to the end of the mitred offset corner.
                When two line segments meet at a sharp angle, a miter join will extend far beyond the original geometry. (and in the extreme case will be infinitely far.) To prevent unreasonable geometry, the mitre limit allows controlling the maximum length of the join corner.
                Corners with a ratio which exceed the limit will be beveled(Default=10)

        """
        self.logger = Logger(feedback)
        self.processor = (
            processor
            if processor
            else AlignerGeometryProcessor(ProcessorConfig(), feedback)
        )
        self.correction_distance = correction_distance
        self.mitre_limit = mitre_limit
        self.max_workers = max_workers

        # PROCESSING DEFAULTS
        # thematic
        # name of the identifier-field of the thematic data (id has to be unique)
        self.name_thematic_id = ID_THEME_FIELD_NAME
        # dictionary to store all thematic geometries to handle
        self.dict_thematic: dict[ThematicId, BaseGeometry] = {}
        # dictionary to store properties of the reference-features (optional)
        self.dict_thematic_properties: dict[ThematicId, dict] = {}
        # Dict to store source-information of the thematic dictionary
        self.dict_thematic_source: dict[ThematicId, str] = {}
        # dictionary to store all unioned thematic geometries
        self.thematic_union = None

        # reference

        # name of the identifier-field of the reference data (id has to be unique,f.e
        # CAPAKEY for GRB-parcels)
        self.name_reference_id = ID_REFERENCE_FIELD_NAME
        # dictionary to store all reference geometries
        self.dict_reference: dict[ThematicId, BaseGeometry] = {}
        # dictionary to store properties of the reference-features (optional)
        self.dict_reference_properties: dict[ThematicId, dict] = {}
        # Dict to store source-information of the reference dictionary
        self.dict_reference_source: dict[ThematicId, str] = {}
        # to save a unioned geometry of all reference polygons; needed for calculation
        # in most OD-strategies
        self.reference_union = None

        # to save a the reference_elements (points and lines) that form the reference borders; needed for networkx-calculation
        # in most OD-strategies
        self.reference_elements = None

        # Coordinate reference system
        # thematic geometries and reference geometries are assumed to be in the same CRS
        # before loading into the Aligner. No CRS-transformation will be performed.
        # When loading data, CRS is expected to be a projected CRS with units in 'meter
        # (m)'.
        # Default EPSG:31370 (Lambert72), alternative: EPSG:3812 (Lambert2008)
        self.CRS = to_crs(crs)
        # TODO The crs that is defined on the aligner is the CRS we are working with. So we expect that the loaded thematic data is in this CRS and also the reference data is in this CRS. Or will be downloaded and transformed to this CRS.
        # At this moment  we expect the units in meter, because all calculation and parameters are based that unit is 'meter'

        # this parameter is used to treat multipolygon as single polygons. So polygons
        # with ID splitter are separately evaluated and merged on result.
        self.multi_as_single_modus = multi_as_single_modus
        self.preserve_topology = preserve_topology
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

    def process(
        self,
        relevant_distances: Iterable[float] = None,
        *,
        dict_thematic_to_process=None,
        max_workers: int = None,
    ) -> AlignerResult:
        """
        Calculates the resulting dictionaries for thematic data based on a series of
            relevant distances.

        Args:
            dict_thematic_to_process: the dictionary with the thematic geometries to 'predict'. Default is None, so all thematic geometries inside the aligner will be processed.
            relevant_distances (Iterable[float]): A series of relevant distances
                (in meters) to process
            max_workers (int, optional): Amount of workers that is used in ThreadPoolExecutor (for parallel execution) when processing objects for multiple relevant distances. (default None). If set to -1, no parallel exececution is used.

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
            raise ValueError("provide at least 1 relevant distance")
        if dict_thematic_to_process is None:
            dict_thematic_to_process = self.dict_thematic
        if max_workers is None:
            max_workers = self.max_workers

        self.logger.feedback_debug("Process series" + str(relevant_distances))

        dict_series: dict[ThematicId, dict[float, ProcessResult | None]] = {}
        futures = {}

        dict_multi_as_single = {}
        topo_thematic = None
        dict_thematic_topo_geoms = None

        if self.multi_as_single_modus:
            dict_thematic_to_process, dict_multi_as_single = multi_to_singles(
                dict_thematic_to_process
            )

        if self.preserve_topology:
            # self.max_workers =-1
            # self.logger.feedback_info("max_workers set to -1 when using 'preserve_topology'")
            dict_thematic_to_process, topo_thematic, dict_thematic_topo_geoms = (
                generate_topo(dict_thematic_to_process)
            )

        def process_geom_for_rd(geometry, relevant_distance):
            return self.processor.process(
                correction_distance=self.correction_distance,
                dict_reference=self.dict_reference,
                mitre_limit=self.mitre_limit,
                reference_elements=self._get_reference_elements(),
                reference_items=self.reference_items,
                reference_tree=self.reference_tree,
                reference_union=self._get_reference_union(),
                input_geometry=geometry,
                relevant_distance=relevant_distance,
            )

        def run_process(executor: ThreadPoolExecutor = None):
            for thematic_id, geom in dict_thematic_to_process.items():
                self.logger.feedback_info(
                    f"thematic id {str(thematic_id)} processed with "
                    f"relevant distances (m) [{str(relevant_distances)}]"
                )
                dict_series[thematic_id] = {}
                for rd in relevant_distances:
                    try:
                        fn = process_geom_for_rd
                        if executor:
                            futures[(thematic_id, rd)] = executor.submit(fn, geom, rd)
                        else:
                            dict_series[thematic_id][rd] = fn(geom, rd)
                    except ValueError as e:
                        self.logger.feedback_warning(
                            f"error for thematic id {str(thematic_id)} processed with "
                            f"relevant distances (m) [{str(relevant_distances)}]"
                        )
                        dict_series[thematic_id][rd] = None
                        self.logger.feedback_warning(str(e))

        if max_workers == -1:
            run_process()
        else:
            with ThreadPoolExecutor(max_workers) as executor:
                run_process(executor)
                self.logger.feedback_debug("waiting all started RD calculations")
                wait(list(futures.values()))
                for (key, rd), future in futures.items():
                    dict_series[key][rd] = future.result()

        if self.preserve_topology:
            dict_series = dissolve_topo(
                dict_series,
                dict_thematic_topo_geoms,
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
                original_geometry_length = len(original_geometry.geoms)  # noqa
            except:
                original_geometry_length = 1
            for relevant_distance, process_result in dict_dist_results.items():
                resulting_geom = process_result["result"]
                try:
                    resulting_geometry_length = len(resulting_geom.geoms)  # noqa
                except:
                    resulting_geometry_length = 1
                if original_geometry_length != resulting_geometry_length:
                    msg = "Difference in amount of geometries"
                    self.logger.feedback_debug(msg)
                    process_result["properties"][REMARK_FIELD_NAME] = (
                        process_result["properties"][REMARK_FIELD_NAME] + " | " + msg
                    )

        self.logger.feedback_info(
            "End of processing series: " + str(relevant_distances)
        )

        return AlignerResult(dict_series)

    def predict(
        self,
        relevant_distances=None,
        *,
        dict_thematic=None,
        diff_metric=None,
    ) -> AlignerResult:
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
            relevant_distances (np.ndarray, optional): A NumPy array of distances to
              be analyzed. Defaults to np.arange(0.1, 5.05, 0.1).
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
        if relevant_distances is None:
            relevant_distances = [
                round(k, RELEVANT_DISTANCE_DECIMALS)
                for k in np.arange(0, 310, 10, dtype=int) / 100
            ]
        rd_prediction = list(relevant_distances)
        max_relevant_distance = max(rd_prediction)
        cvg_ratio = coverage_ratio(values=relevant_distances, min_val=0, bin_count=10)
        cvg_ratio_threshold = 0.75
        # cvg_ratio: indication of the rd values can be used to make a brdr_prediction_score. When there is enough coverage of predictions to determine a prediction_score we also add 0 and a long-range value(+1).
        # Otherwise we only add a short-range value (+0.1) to check for stability
        if cvg_ratio > cvg_ratio_threshold:
            rd_prediction.append(round(0, RELEVANT_DISTANCE_DECIMALS))
            rd_prediction.append(
                round(max_relevant_distance + 0.1, RELEVANT_DISTANCE_DECIMALS)
            )
            rd_prediction.append(
                round(max_relevant_distance + 1, RELEVANT_DISTANCE_DECIMALS)
            )
        else:
            rd_prediction.append(
                round(max_relevant_distance + 0.1, RELEVANT_DISTANCE_DECIMALS)
            )
        rd_prediction = list(set(rd_prediction))
        rd_prediction = sorted(rd_prediction)

        process_result = self.process(
            dict_thematic_to_process=dict_thematic,
            relevant_distances=rd_prediction,
        )
        if diff_metric is None:
            diff_metric = self.diff_metric
        diffs_dict = {}
        for theme_id, dict_processresult in process_result.results.items():
            diffs = diffs_from_dict_processresult(
                dict_processresult,
                dict_thematic[theme_id],
                self._get_reference_union(),
                diff_metric=diff_metric,
            )
            diffs_dict[theme_id] = diffs
            if len(diffs) != len(rd_prediction):
                self.logger.feedback_warning(
                    f"Number of computed diffs for thematic element {theme_id} does "
                    f"not match the number of relevant distances."
                )
                continue
            diff_values = list(diffs.values())
            dict_stability = determine_stability(rd_prediction, diff_values)
            for rd in rd_prediction:
                if rd not in relevant_distances:
                    del process_result.results[theme_id][rd]
                    continue
                process_result.results[theme_id][rd]["properties"][STABILITY] = (
                    dict_stability[rd][STABILITY]
                )
                if dict_stability[rd][ZERO_STREAK] is not None:
                    if cvg_ratio > cvg_ratio_threshold:
                        process_result.results[theme_id][rd]["properties"][
                            PREDICTION_SCORE
                        ] = dict_stability[rd][ZERO_STREAK][3]

        self.diffs_dict=diffs_dict
        self.count_predictions(process_result.results)
        return process_result

    def count_predictions(self, dict_predictions: dict[ThematicId, dict[float, ProcessResult]]):
        """
        # Check if the predicted geometries are unique (and remove duplicated predictions)
        """
        #TODO moet deze functie niets teruggeven?
        dict_predictions_unique = defaultdict(dict)
        for theme_id, dist_results in dict_predictions.items():
            dict_predictions_unique[theme_id] = {}
            predicted_geoms_for_theme_id = []
            for rel_dist, processresult in dist_results.items():
                predicted_geom = processresult["result"]
                if not equal_geom_in_array(
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
                    dict_predictions_unique[theme_id][rel_dist] = processresult
                    predicted_geoms_for_theme_id.append(processresult["result"])
                else:
                    self.logger.feedback_info(
                        f"Duplicate prediction found for key {theme_id} at distance {rel_dist}"
                    )
            for dist, processresult in dist_results.items():
                processresult["properties"][PREDICTION_COUNT] = len(predicted_geoms_for_theme_id)

    def evaluate(
        self,
        relevant_distances = None,
        *,
        ids_to_evaluate=None,
        base_formula_field=FORMULA_FIELD_NAME,
        full_strategy=FullStrategy.NO_FULL,
        max_predictions=-1,
        multi_to_best_prediction=True,
    ):
        """
        Compares and evaluate input-geometries (with formula). Attributes are added to evaluate and decide if new
        proposals can be used
        ids_to_evaluate: list with all IDs to evaluate. all other IDs will be unchanged. If None (default), all self.dict_thematic will be evaluated.
        base_formula_field: name of the field where the base_formula is found in the data
        max_predictions: integer that indicates how mstr|int predictions are maximally returned. (-1 indicates all predictions are returned)
        relevant_distances: relevant distances to evaluate
        full_strategy: enum, decided which predictions are kept or prefered based on full-ness of the prediction
        multi_to_best_prediction (default True): Only usable in combination with max_predictions=1. If True (and max_predictions=1), the prediction with highest score will be taken.If False, the original geometry is returned.
        """

        if ids_to_evaluate is None:
            ids_to_evaluate = list(self.dict_thematic.keys())
        if any(id_to_evaluate not in self.dict_thematic.keys() for id_to_evaluate in ids_to_evaluate):
            raise ValueError("not all ids_to_evaluate are found in the thematic data")
        if relevant_distances is None:
            relevant_distances = [
                round(k, RELEVANT_DISTANCE_DECIMALS)
                for k in np.arange(0, 310, 10, dtype=int) / 100
            ]
        dict_affected = {}
        dict_unaffected = {}
        for id_theme, geom in self.dict_thematic.items():
            if id_theme in ids_to_evaluate:
                dict_affected[id_theme] = geom
            else:
                dict_unaffected[id_theme] = geom

        # Features are split up in 2 dicts: affected and unaffected (no_change)
        # The affected features will be split up:
        #   *No prediction available
        #   *Predictions available

        # AFFECTED
        prediction_result = self.predict(
            dict_thematic=dict_affected,
            relevant_distances=relevant_distances,
            diff_metric=self.diff_metric,
        )
        dict_affected_predictions = prediction_result.get_results(AlignerResultType.PREDICTIONS)
        dict_predictions_evaluated = {}

        for theme_id in dict_affected.keys():
            dict_predictions_evaluated[theme_id] = {}
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
                props[REMARK_FIELD_NAME] = "no predictions available, original returned"
                dict_predictions_evaluated[theme_id][relevant_distance] = {
                    "result": dict_affected[theme_id],
                    "properties": props,
                }
                continue

            # When there are predictions available
            dict_predictions_results = dict_affected_predictions[theme_id]
            scores = []
            distances = []
            predictions = []
            formula_match = False
            for dist in sorted(dict_predictions_results.keys()):
                props = self._evaluate(
                    id_theme=theme_id,
                    geom_predicted=dict_predictions_results[dist]["result"],
                    base_formula_field=base_formula_field,
                )
                props.update(dict_affected_predictions[theme_id][dist]["properties"])

                full = props[FULL_ACTUAL_FIELD_NAME]
                if full_strategy == FullStrategy.ONLY_FULL and not full:
                    # this prediction is ignored
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
                    # this prediction has a equality based on formula so the rest is not checked anymore
                    formula_match = True
                    props[PREDICTION_SCORE] = 100
                    scores = []
                    distances = []
                    predictions = []
                    scores.append(props[PREDICTION_SCORE])
                    distances.append(dist)
                    dict_affected_predictions[theme_id][dist]["properties"] = props
                    predictions.append(dict_affected_predictions[theme_id][dist])
                    continue
                if full:
                    if full_strategy != FullStrategy.NO_FULL:
                        props[EVALUATION_FIELD_NAME] = (
                            Evaluation.TO_CHECK_PREDICTION_FULL
                        )
                        prediction_score = props[PREDICTION_SCORE] + 50
                        if prediction_score > 100:
                            prediction_score = 100
                        props[PREDICTION_SCORE] = prediction_score
                    else:
                        props[EVALUATION_FIELD_NAME] = (
                            Evaluation.TO_CHECK_PREDICTION_MULTI_FULL
                        )

                scores.append(props[PREDICTION_SCORE])
                distances.append(dist)
                dict_affected_predictions[theme_id][dist]["properties"] = props
                predictions.append(dict_affected_predictions[theme_id][dist])

            # get max amount of best-scoring predictions
            best_ix = sorted(range(len(scores)), reverse=True, key=lambda i: scores[i])
            len_best_ix = len(best_ix)

            if not formula_match:
                # if there is only one prediction left,  evaluation is set to PREDICTION_UNIQUE_FULL
                if len_best_ix == 1 and not formula_match:
                    props = predictions[0]["properties"]
                    if (
                        FULL_ACTUAL_FIELD_NAME in props
                        and props[FULL_ACTUAL_FIELD_NAME]
                    ):
                        predictions[0]["properties"][
                            EVALUATION_FIELD_NAME
                        ] = Evaluation.PREDICTION_UNIQUE_FULL
                    else:
                        predictions[0]["properties"][
                            EVALUATION_FIELD_NAME
                        ] = Evaluation.PREDICTION_UNIQUE

                # if there are multiple predictions, but we want only one and we ask for the original
                if (
                    len_best_ix > 1
                    and max_predictions == 1
                    and not multi_to_best_prediction
                    and not formula_match
                ):
                    relevant_distance = round(0, RELEVANT_DISTANCE_DECIMALS)
                    props[EVALUATION_FIELD_NAME] = Evaluation.TO_CHECK_ORIGINAL
                    props[PREDICTION_SCORE] = -1
                    props[REMARK_FIELD_NAME] = "multiple predictions, original returned"
                    dict_predictions_evaluated[theme_id][relevant_distance] = {
                        "result": dict_affected[theme_id],
                        "properties": props,
                    }
                    continue

            if max_predictions > 0 and len_best_ix > max_predictions:
                best_ix = best_ix[:max_predictions]
            if len(best_ix) > 0:
                for ix in best_ix:
                    distance = distances[ix]
                    prediction = predictions[ix]
                    dict_predictions_evaluated[theme_id][distance] = prediction
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
                props[REMARK_FIELD_NAME] = "no prediction, original returned"
                dict_predictions_evaluated[theme_id][relevant_distance] = {
                    "result": dict_affected[theme_id],
                    "properties": props,
                }

        # UNAFFECTED
        relevant_distance = round(0, RELEVANT_DISTANCE_DECIMALS)
        for theme_id, geom in dict_unaffected.items():
            dict_predictions_evaluated[theme_id] = {}
            props = self._evaluate(
                id_theme=theme_id,
                geom_predicted=geom,
                base_formula_field=base_formula_field,
            )
            props[EVALUATION_FIELD_NAME] = Evaluation.NO_CHANGE
            props[PREDICTION_SCORE] = -1
            props[REMARK_FIELD_NAME] = (
                "Unaffected (no change) --> original geometry returned"
            )
            dict_predictions_evaluated[theme_id][relevant_distance] = {
                "result": geom,
                "properties": props,
            }

        return AlignerResult(dict_predictions_evaluated)

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
        dict_series (dict): A dictionary where keys are thematic IDs and values are dictionaries mapping relevant distances to ProcessResult objects.
        dict_thematic (dict): A dictionary where keys are thematic IDs and values are BaseGeometry objects representing the original geometries.
        diff_metric (DiffMetric, optional): The metric to use for calculating differences. Default is DiffMetric.CHANGES_AREA.

        Returns:
        dict: A dictionary where keys are thematic IDs and values are dictionaries mapping relevant distances to calculated difference metrics.
        """
        if dict_processresults is None:
            raise ValueError("dict_processresults is required")
        if dict_thematic is None:
            dict_thematic = self.dict_thematic
        diffs = {}
        for key in dict_thematic:
            diffs[key] = diffs_from_dict_processresult(
                dict_processresult=dict_processresults[key],
                geom_thematic=dict_thematic[key],
                reference_union=self._get_reference_union(),
                diff_metric=diff_metric,
            )

        return diffs

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

    def _get_reference_union(self) -> BaseGeometry:
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
