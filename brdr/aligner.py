import hashlib
import inspect
import json
import os
import uuid
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import wait
from copy import deepcopy
from datetime import datetime
from typing import Iterable, Dict, Any, Optional

import numpy as np
from shapely import make_valid
from shapely import to_geojson
from shapely.geometry.base import BaseGeometry

from brdr import __version__
from brdr.configs import ProcessorConfig
from brdr.constants import (
    DATE_FORMAT,
    SYMMETRICAL_AREA_CHANGE,
    SYMMETRICAL_AREA_PERCENTAGE_CHANGE,
    AREA_CHANGE,
    AREA_PERCENTAGE_CHANGE,
    LENGTH_CHANGE,
    LENGTH_PERCENTAGE_CHANGE,
)
from brdr.constants import DEFAULT_CRS
from brdr.constants import DIFF_AREA_FIELD_NAME
from brdr.constants import DIFF_METRIC
from brdr.constants import DIFF_PERCENTAGE_FIELD_NAME
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
from brdr.enums import FullReferenceStrategy
from brdr.enums import ProcessRemark
from brdr.feature_data import AlignerFeatureCollection
from brdr.geometry_utils import buffer_neg, geometric_equality
from brdr.geometry_utils import buffer_pos
from brdr.geometry_utils import safe_difference
from brdr.geometry_utils import safe_intersection
from brdr.geometry_utils import safe_unary_union
from brdr.geometry_utils import to_crs
from brdr.loader import Loader
from brdr.logger import Logger
from brdr.processor import AlignerGeometryProcessor
from brdr.processor import BaseProcessor
from brdr.typings import ProcessResult
from brdr.typings import ThematicId
from brdr.utils import (
    coverage_ratio,
    recursive_stepwise_interval_check,
    create_full_interpolated_dataset,
)
from brdr.utils import unary_union_result_dict
from brdr.utils import determine_stability
from brdr.utils import geojson_from_dict
from brdr.utils import get_geojsons_from_process_results
from brdr.utils import get_geometry_difference_metrics_from_processresult
from brdr.utils import get_geometry_difference_metrics_from_processresults
from brdr.utils import is_brdr_formula
from brdr.utils import write_geojson


###################
class AlignerResult:
    """
    Stores the results from the alignment process for different themes and relevant distances.

    This class is responsible for storing, enriching, and filtering the
    processed geometries and their associated metrics.

    Attributes
    ----------
    metadata : Any | None
        Optional metadata about the overall alignment operation.
    results : dict[ThematicId, dict[float, ProcessResult | None]]
        The raw processing results, structured by theme ID and relevant distance.
    """

    metadata: Optional[Any] = None

    def __init__(
        self,
        process_results: dict[ThematicId, dict[float, ProcessResult | None]],
    ):
        """
        Initializes an AlignerResult instance.

        Parameters
        ----------
        process_results : dict[ThematicId, dict[float, ProcessResult | None]]
            A nested dictionary containing the raw processing results.
            The outer key is the theme ID, the inner key is the 'relevant distance'
            used for the calculation.
        """
        self.results = process_results

    def get_results(
        self,
        aligner: "Aligner",
        result_type: AlignerResultType = AlignerResultType.PROCESSRESULTS,
    ) -> dict[ThematicId, dict[float, ProcessResult | None]]:
        """
        Retrieves the processing results, enriching them with geometric difference metrics,
        and filtering them based on the requested result type.

        This method mutates the 'properties' of the stored ProcessResult objects by
        adding several calculated geometric difference metrics.

        Parameters
        ----------
        aligner : Aligner
            The 'Aligner' object containing the original geometries and comparison logic.
        result_type : AlignerResultType, optional
            The type of results to return. Options are:
            - PROCESSRESULTS: All enriched results (default).
            - PREDICTIONS: Only results containing a prediction score.
            - EVALUATED_PREDICTIONS: Only results that have been evaluated.

        Returns
        -------
        dict[ThematicId, dict[float, ProcessResult | None]]
            The filtered and enriched dictionary of processing results.

        Raises
        ------
        ValueError
            If an unknown `result_type` is provided.
        """
        for theme_id, results_dict in self.results.items():
            original_geometry = aligner.dict_thematic[theme_id]
            nr_calculations = len(results_dict)
            try:
                original_geometry_length = len(original_geometry.geoms)
            except:
                original_geometry_length = 1

            for relevant_distance, process_result in results_dict.items():
                if process_result is None:
                    continue

                properties = process_result["properties"]
                # Adding extra properties
                properties[ID_THEME_FIELD_NAME] = theme_id
                properties[NR_CALCULATION_FIELD_NAME] = nr_calculations
                properties[RELEVANT_DISTANCE_FIELD_NAME] = relevant_distance

                # --- Adding Geometric Difference Metrics ---
                properties[SYMMETRICAL_AREA_CHANGE] = (
                    get_geometry_difference_metrics_from_processresult(
                        process_result,
                        aligner.dict_thematic[theme_id],
                        None,
                        DIFF_METRIC.SYMMETRICAL_AREA_CHANGE,
                    )
                )
                properties[SYMMETRICAL_AREA_PERCENTAGE_CHANGE] = (
                    get_geometry_difference_metrics_from_processresult(
                        process_result,
                        aligner.dict_thematic[theme_id],
                        None,
                        DIFF_METRIC.SYMMETRICAL_AREA_PERCENTAGE_CHANGE,
                    )
                )
                properties[AREA_CHANGE] = (
                    get_geometry_difference_metrics_from_processresult(
                        process_result,
                        aligner.dict_thematic[theme_id],
                        None,
                        DIFF_METRIC.AREA_CHANGE,
                    )
                )
                properties[AREA_PERCENTAGE_CHANGE] = (
                    get_geometry_difference_metrics_from_processresult(
                        process_result,
                        aligner.dict_thematic[theme_id],
                        None,
                        DIFF_METRIC.AREA_PERCENTAGE_CHANGE,
                    )
                )
                properties[LENGTH_CHANGE] = (
                    get_geometry_difference_metrics_from_processresult(
                        process_result,
                        aligner.dict_thematic[theme_id],
                        None,
                        DIFF_METRIC.LENGTH_CHANGE,
                    )
                )
                properties[LENGTH_PERCENTAGE_CHANGE] = (
                    get_geometry_difference_metrics_from_processresult(
                        process_result,
                        aligner.dict_thematic[theme_id],
                        None,
                        DIFF_METRIC.LENGTH_PERCENTAGE_CHANGE,
                    )
                )

                # --- Checking for Geometry Count Change ---
                resulting_geom = process_result["result"]
                try:
                    resulting_geometry_length = len(resulting_geom.geoms)
                except:
                    resulting_geometry_length = 1

                if original_geometry_length != resulting_geometry_length:
                    remark = ProcessRemark.CHANGED_AMOUNT_GEOMETRIES
                    remarks = properties.get(REMARK_FIELD_NAME, [])
                    remarks.append(remark)
                    properties[REMARK_FIELD_NAME] = remarks

        if result_type == AlignerResultType.PROCESSRESULTS:
            return self.results
        elif result_type == AlignerResultType.PREDICTIONS:
            # return all theme ids with only the relevant distances where the ProcessResult has property PREDICTION_SCORE:
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
            # return all theme ids with only the relevant distances where the ProcessResult has property EVALUATION_FIELD_NAME:
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
        aligner: "Aligner",
        result_type: AlignerResultType = AlignerResultType.PROCESSRESULTS,
        formula: bool = False,
        attributes: bool = False,
    ) -> Dict[str, Any]:
        """
        Returns the processing results as a GeoJSON FeatureCollection.

        The GeoJSON contains the resulting geometries for all 'serial' relevant distances
        within the selected result type, optionally enriched with original attributes and comparison formulas.

        Parameters
        ----------
        aligner : Aligner
            The 'Aligner' object, necessary to retrieve the CRS, ID fields, and comparison formulas.
        result_type : AlignerResultType, optional
            The type of results to export (default is PROCESSRESULTS).
        formula : bool, optional
            If True, the descriptive comparison formula is added as an attribute to the result.
            Defaults to False.
        attributes : bool, optional
            If True, the original attributes/properties of the thematic objects are
            added to the result. Defaults to False.

        Returns
        -------
        dict
            A dictionary representing a GeoJSON FeatureCollection.

        Raises
        ------
        ValueError
            If the results are empty (`self.results` is None or empty).
        """

        if self.results is None or self.results == {}:
            raise ValueError("Empty results: No calculated results to export.")

        results = self.get_results(aligner=aligner, result_type=result_type)
        prop_dictionary: defaultdict[Any, dict[Any, dict]] = defaultdict(dict)

        for theme_id, results_dict in results.items():
            prop_dictionary[theme_id] = {}
            for relevant_distance, process_result in results_dict.items():
                if process_result is None:
                    continue

                prop_dictionary[theme_id][relevant_distance] = {}

                # Adding original attributes
                if attributes and theme_id in aligner.dict_thematic_properties.keys():
                    for attr, value in aligner.dict_thematic_properties[
                        theme_id
                    ].items():
                        prop_dictionary[theme_id][relevant_distance][attr] = value

                # Adding formula
                if formula:
                    result = process_result["result"]
                    formula_result = aligner.compare_to_reference(result)
                    prop_dictionary[theme_id][relevant_distance][FORMULA_FIELD_NAME] = (
                        json.dumps(formula_result)
                    )

        return get_geojsons_from_process_results(
            results,
            crs=aligner.CRS,
            id_field=aligner.name_thematic_id,
            series_prop_dict=prop_dictionary,
        )

    def save_results(
        self,
        aligner,
        path,
        result_type: AlignerResultType = AlignerResultType.PROCESSRESULTS,
        formula=False,
        attributes=False,
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
            aligner=aligner,
            result_type=result_type,
            formula=formula,
            attributes=attributes,
        )
        for name, fc in fcs.items():
            write_geojson(
                os.path.join(path, result_type.value + "_" + name + ".geojson"), fc
            )


def aligner_metadata_decorator(f):
    def inner_func(aligner, *args, **kwargs):
        assert isinstance(aligner, Aligner)
        response = f(aligner, *args, **kwargs)
        if aligner.log_metadata:
            # generate uuid for actuation
            actuation_id = "brdrid:actuations/" + uuid.uuid4().hex
            processor_id = aligner.processor.processor_id.value
            processor_name = type(aligner.processor).__name__
            reference_data = aligner.reference_data
            reference_features = reference_data.features.values()
            reference_geometries = [
                {
                    "id": feature.brdr_id,
                    "type": feature.geometry.geom_type,
                    "version_date": reference_data.source.get(VERSION_DATE, ""),
                }
                for feature in reference_features
            ]
            stack = inspect.stack()
            f_locals = stack[0][0].f_locals

            for thematic_id, rd, result in [
                (thematic_id, rd, res)
                for thematic_id, rd_res in response.results.items()
                for rd, res in rd_res.items()
            ]:
                thematic_feature = aligner.thematic_data.features[thematic_id]
                feature_of_interest_id = thematic_feature.brdr_id  # TODO
                result_hash = hashlib.sha256(result["result"].wkt.encode()).hexdigest()
                result["metadata"] = {
                    "id": actuation_id,
                    "type": "sosa:Actuation",
                    "reference_geometries": reference_geometries,
                    "changes": "geo:hasGeometry",
                    "sosa:hasFeatureOfInterest": {"id": feature_of_interest_id},
                    "result": f"brdrid:geoms/{result_hash}",  # id is sha265 hash of wkt
                    "procedure": {
                        "id": processor_id,
                        "implementedBy": processor_name,
                        "type": "sosa:Procedure",
                        "ssn:hasInput": [
                            {
                                "id": "brdr:relevante_afstand",
                                "type": "ssn:Input",
                                "input_value": {"type": "xsd:integer", "value": rd},
                            },
                        ],
                    },
                }

        return response

    return inner_func

    # TODO what about the Aligner-parameters; AlignerConfig-class?
class Aligner:
    """
    Compares and aligns thematic geospatial data against a set of reference data.

    The Aligner manages the loading of both thematic and reference data, configures
    the geometric processing rules, and executes the alignment, prediction, and
    evaluation logic across a series of relevant distances.

    Attributes
    ----------
    logger : Logger
        Instance for logging feedback and information.
    log_metadata : bool
        If True, metadata about the actuation is logged in the results.
    processor : BaseProcessor | AlignerGeometryProcessor
        The geometric processor used for alignment calculations.
    correction_distance : float
        Distance used in buffer operations to remove slivers (technical correction).
    mitre_limit : int
        Parameter for the buffer operation to control the maximum length of join corners.
    max_workers : int | None
        The maximum number of workers for parallel execution (ThreadPoolExecutor).
    CRS : str
        The Coordinate Reference System (CRS) being used (e.g., 'EPSG:31370').
    name_thematic_id : str
        Name of the identifier field for thematic data.
    dict_thematic : dict[ThematicId, BaseGeometry]
        Dictionary storing all thematic geometries.
    dict_thematic_properties : dict[ThematicId, dict]
        Dictionary storing properties of thematic features.
    diff_metric : DiffMetric
        The metric used to measure differences between geometries (e.g., area change).
    reference_data : AlignerFeatureCollection | None
        Loaded collection of reference features.
    thematic_data : AlignerFeatureCollection | None
        Loaded collection of thematic features.
    """

    def __init__(
        self,
        *,
        feedback=None,
        processor: BaseProcessor = None,
        crs=DEFAULT_CRS,
        correction_distance=0.01,
        diff_metric=DIFF_METRIC,
        mitre_limit=10,
        max_workers=None,
        log_metadata=True,
    ):
        """
        Initializes the Aligner object.

        Parameters
        ----------
        feedback : object, optional
            Feedback object used for logging (e.g., in a QGIS environment). Defaults to None.
        processor : BaseProcessor, optional
            The geometric processor instance. If None, AlignerGeometryProcessor is used.
        crs : str, optional
            Coordinate Reference System (CRS) of the data (default is EPSG:31370).
            Expected to be a projected CRS with units in 'meter (m)'.
        correction_distance : float, optional
            Distance (in meters) used in a positive/negative buffer operation to remove slivers
            or define geometric equality. Default is 0.01m (1cm).
        diff_metric : DiffMetric, optional
            Metric used to measure differences in alignment/prediction (e.g., symmetrical area change).
            This value may be overwritten internally for point/line geometries.
        mitre_limit : int, optional
            Mitre ratio limit used in buffer operations. Corners exceeding this ratio
            will be beveled to prevent excessive extensions. Defaults to 10.
        max_workers : int | None, optional
            Maximum number of workers for parallel processing. If None, default system
            concurrency is used. If -1, no parallel execution is used.
        log_metadata : bool, optional
            If True, actuation metadata is generated and attached to results (following
            the SOSA/SSN standard). Defaults to True.
        """
        self.logger = Logger(feedback)
        self.log_metadata = log_metadata
        self.processor = (
            processor
            if processor
            else AlignerGeometryProcessor(ProcessorConfig(), feedback)
        )
        self.correction_distance = correction_distance
        self.mitre_limit = mitre_limit
        self.max_workers = max_workers

        # PROCESSING DEFAULTS (Internal state variables)
        self.name_thematic_id = ID_THEME_FIELD_NAME
        self.dict_thematic: dict[ThematicId, BaseGeometry] = {}
        self.dict_thematic_properties: dict[ThematicId, dict] = {}
        self.dict_thematic_source: dict[ThematicId, str] = {}

        self.name_reference_id = ID_REFERENCE_FIELD_NAME
        # The CRS is the working CRS for all calculations (assumed to be projected/in meters)
        self.CRS = to_crs(crs)

        self.diff_metric = diff_metric
        self.logger.feedback_info("Aligner initialized")

        self.reference_data: AlignerFeatureCollection | None = None
        self.thematic_data: AlignerFeatureCollection | None = None

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
        self.thematic_data = loader.load_data_as_feature_collection()


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
        self.reference_data = loader.load_data_as_feature_collection()
        self.reference_data.is_reference = True

    @aligner_metadata_decorator
    def process(
        self,
        relevant_distances: Iterable[float] = None,
        *,
        dict_thematic=None,
        max_workers: int = None,
    ) -> AlignerResult:
        """
        Calculates the resulting dictionaries for thematic data based on a series of
            relevant distances.

        Args:
            dict_thematic: the dictionary with the thematic geometries to 'predict'. Default is None, so all thematic geometries inside the aligner will be processed.
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
        if dict_thematic is None:
            dict_thematic = self.dict_thematic
        if any(
            id_to_process not in self.dict_thematic.keys()
            for id_to_process in dict_thematic.keys()
        ):
            raise ValueError("not all ids are found in the thematic data")
        if max_workers is None:
            max_workers = self.max_workers

        self.logger.feedback_debug("Process series" + str(relevant_distances))

        process_results: dict[ThematicId, dict[float, ProcessResult | None]] = {}
        futures = {}

        def process_geom_for_rd(geometry, relevant_distance):
            return self.processor.process(
                correction_distance=self.correction_distance,
                reference_data=self.reference_data,
                mitre_limit=self.mitre_limit,
                input_geometry=geometry,
                relevant_distance=relevant_distance,
                dict_thematic=self.dict_thematic,
            )

        def run_process(executor: ThreadPoolExecutor = None):
            for thematic_id, geom in dict_thematic.items():
                self.logger.feedback_info(
                    f"thematic id {str(thematic_id)} processed with "
                    f"relevant distances (m) [{str(relevant_distances)}]"
                )
                process_results[thematic_id] = {}
                for rd in relevant_distances:
                    try:
                        fn = process_geom_for_rd
                        if executor:
                            futures[(thematic_id, rd)] = executor.submit(fn, geom, rd)
                        else:
                            process_results[thematic_id][rd] = fn(geom, rd)
                    except ValueError as e:
                        self.logger.feedback_warning(
                            f"error for thematic id {str(thematic_id)} processed with "
                            f"relevant distances (m) [{str(relevant_distances)}]"
                        )
                        process_results[thematic_id][rd] = None
                        self.logger.feedback_warning(str(e))

        if max_workers == -1:
            run_process()
        else:
            with ThreadPoolExecutor(max_workers) as executor:
                run_process(executor)
                self.logger.feedback_debug("waiting all started RD calculations")
                wait(list(futures.values()))
                for (key, rd), future in futures.items():
                    process_results[key][rd] = future.result()

        self.logger.feedback_info(
            "End of processing series: " + str(relevant_distances)
        )

        return AlignerResult(process_results)

    def predict(
        self,
        relevant_distances=None,
        *,
        dict_thematic=None,
        diff_metric=None,
        process_all_at_once=True,

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
            process_all_at_once=True
                #True: All calculations are done for all relevant distances. This seems to be faster than the other method of doing calculations for some intermediate points and copy result when the same.
                #False: Uses the logic to only calculate intermediate points and copy if equal, otherwise extra points are calculated. Until the full range is filled.

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
        if dict_thematic is None:  # or dict_thematic =={}
            dict_thematic = self.dict_thematic
        if any(
            id_to_predict not in self.dict_thematic.keys()
            for id_to_predict in dict_thematic.keys()
        ):
            raise ValueError("not all ids are found in the thematic data")
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
        if process_all_at_once:
            aligner_result = self.process(
                dict_thematic=dict_thematic,
                relevant_distances=rd_prediction,
            )
            process_results =  aligner_result.results

        else:
            process_results = {}
            futures={}

            def _predict_geom(theme_id,geom,_rd_prediction):

                def _check_interval_stability(
                    res_start: ProcessResult, res_end: ProcessResult
                ) -> bool:

                    return geometric_equality(
                        res_start["result"],
                        res_end["result"],
                        correction_distance=self.correction_distance,
                        mitre_limit=self.mitre_limit,
                    )

                def _process_result(_theme_id, _relevant_distances):
                    aligner_result = self.process(
                        dict_thematic={_theme_id: geom},
                        relevant_distances=_relevant_distances,
                        max_workers=self.max_workers)
                    return aligner_result.results[_theme_id][_relevant_distances[0]]

                def _process_wrapper(x: float):
                    return _process_result(_theme_id=theme_id, _relevant_distances=[x])

                non_stable_points, cache = recursive_stepwise_interval_check(
                    f_value=_process_wrapper,
                    f_condition=_check_interval_stability,
                    discrete_list=_rd_prediction,
                    initial_sample_size=3,
                )
                interpolated_cache = create_full_interpolated_dataset(_rd_prediction, cache)
                return interpolated_cache

            def run_prediction(executor: ThreadPoolExecutor = None):

                for theme_key, geom in dict_thematic.items():

                    self.logger.feedback_info(
                        f"thematic id {str(theme_key)} predictions calculated with "
                        f"relevant distances (m) [{str(relevant_distances)}]"
                    )
                    process_results[theme_key] = {}
                    try:
                        fn = _predict_geom
                        if executor:
                            futures[theme_key] = executor.submit(fn, theme_key,geom,rd_prediction)
                        else:
                            process_results[theme_key] = fn(theme_key,geom,rd_prediction)
                    except ValueError as e:
                        self.logger.feedback_warning(
                            f"error for thematic id {str(theme_key)} processed with "
                            f"relevant distances (m) [{str(relevant_distances)}]"
                        )
                        process_results[theme_key] = None
                        self.logger.feedback_warning(str(e))

            if self.max_workers == -1:
                run_prediction()
            else:
                with ThreadPoolExecutor(self.max_workers) as prediction_executor:
                    run_prediction(prediction_executor)
                    self.logger.feedback_debug("waiting all started thematic calculations")
                    wait(list(futures.values()))
                    for thematic_id, future in futures.items():
                        process_results[thematic_id] = future.result()

        if diff_metric is None:
            diff_metric = self.diff_metric
        diffs_dict = {}
        for theme_id, process_result in process_results.items():
            diffs = get_geometry_difference_metrics_from_processresults(
                process_result,
                dict_thematic[theme_id],
                self.reference_data.union,
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
            prediction_count = 0
            for rd in rd_prediction:
                if rd not in relevant_distances:
                    del process_results[theme_id][rd]
                    continue
                process_results[theme_id][rd]["properties"][STABILITY] = dict_stability[
                    rd
                ][STABILITY]
                if dict_stability[rd][ZERO_STREAK] is not None:
                    if cvg_ratio > cvg_ratio_threshold:
                        prediction_count += 1
                        process_results[theme_id][rd]["properties"][
                            PREDICTION_SCORE
                        ] = dict_stability[rd][ZERO_STREAK][3]
            for rd, process_result in process_results[theme_id].items():
                if PREDICTION_SCORE in process_result["properties"]:
                    process_result["properties"][PREDICTION_COUNT] = prediction_count

        self.diffs_dict = diffs_dict
        return AlignerResult(process_results)

    def evaluate(
        self,
        relevant_distances=None,
        *,
        dict_thematic=None,
        base_formula_field=FORMULA_FIELD_NAME,
        full_reference_strategy=FullReferenceStrategy.NO_FULL_REFERENCE,
        max_predictions=-1,
        multi_to_best_prediction=True,
        process_all_at_once=True,
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

        if dict_thematic is None:
            ids_to_evaluate = list(self.dict_thematic.keys())
        else:
            ids_to_evaluate = list(dict_thematic.keys())
        if any(
            id_to_evaluate not in self.dict_thematic.keys()
            for id_to_evaluate in ids_to_evaluate
        ):
            raise ValueError("not all ids are found in the thematic data")
        if relevant_distances is None:
            relevant_distances = [
                round(k, RELEVANT_DISTANCE_DECIMALS)
                for k in np.arange(0, 310, 10, dtype=int) / 100
            ]
        dict_evaluated = {}
        dict_not_evaluated = {}
        for id_theme, geom in self.dict_thematic.items():
            if id_theme in ids_to_evaluate:
                dict_evaluated[id_theme] = geom
            else:
                dict_not_evaluated[id_theme] = geom

        # Features are split up in 2 dicts: EVALUATED and NOT_EVALUATED (original returned)
        # The evaluated features will be split up:
        #   *No prediction available
        #   *Predictions available

        # PART 1: EVALUATED
        aligner_result = self.predict(
            dict_thematic=dict_evaluated,
            relevant_distances=relevant_distances,
            diff_metric=self.diff_metric,
            process_all_at_once=process_all_at_once,
        )
        process_results_evaluated = aligner_result.get_results(aligner=self)
        process_results_evaluated_predictions = aligner_result.get_results(aligner=self,result_type=AlignerResultType.PREDICTIONS
        )
        process_results_evaluated = deepcopy(process_results_evaluated)

        for theme_id in dict_evaluated.keys():
            if theme_id not in process_results_evaluated_predictions.keys():
                # No predictions available
                relevant_distance = round(0, RELEVANT_DISTANCE_DECIMALS)
                props = self._evaluate(
                    id_theme=theme_id,
                    geom_predicted=dict_evaluated[theme_id],
                    base_formula_field=base_formula_field,
                )
                props[EVALUATION_FIELD_NAME] = Evaluation.TO_CHECK_NO_PREDICTION
                props[PREDICTION_COUNT] = 0
                props[PREDICTION_SCORE] = -1
                if REMARK_FIELD_NAME in props:
                    remarks = props[REMARK_FIELD_NAME]
                else:
                    remarks = []
                remarks.append(ProcessRemark.NO_PREDICTION_ORIGINAL_RETURNED)
                props[REMARK_FIELD_NAME] = remarks
                process_results_evaluated[theme_id][relevant_distance] = unary_union_result_dict({
                    "result": dict_evaluated[theme_id],
                    "properties": props,
                })
                continue

            # When there are predictions available
            dict_predictions_results = process_results_evaluated_predictions[theme_id]
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
                props.update(process_results_evaluated_predictions[theme_id][dist]["properties"])

                full = props[FULL_ACTUAL_FIELD_NAME]
                if (
                    full_reference_strategy == FullReferenceStrategy.ONLY_FULL_REFERENCE
                    and not full
                ):
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
                    process_results_evaluated_predictions[theme_id][dist]["properties"] = props
                    predictions.append(process_results_evaluated_predictions[theme_id][dist])
                    continue
                if full:
                    if full_reference_strategy != FullReferenceStrategy.NO_FULL_REFERENCE:
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
                process_results_evaluated_predictions[theme_id][dist]["properties"] = props
                predictions.append(process_results_evaluated_predictions[theme_id][dist])

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
                    if REMARK_FIELD_NAME in props:
                        remarks = props[REMARK_FIELD_NAME]
                    else:
                        remarks = []
                    remarks.append(ProcessRemark.MULTIPLE_PREDICTIONS_ORIGINAL_RETURNED)
                    props[REMARK_FIELD_NAME] = remarks
                    process_results_evaluated[theme_id][relevant_distance] = unary_union_result_dict({
                        "result": dict_evaluated[theme_id],
                        "properties": props,
                    })
                    continue

            if max_predictions > 0 and len_best_ix > max_predictions:
                best_ix = best_ix[:max_predictions]
            if len(best_ix) > 0:
                for ix in best_ix:
                    distance = distances[ix]
                    prediction = predictions[ix]
                    process_results_evaluated[theme_id][distance] = prediction
            else:
                # #when no evaluated predictions, the original is returned
                relevant_distance = round(0, RELEVANT_DISTANCE_DECIMALS)
                props = self._evaluate(
                    id_theme=theme_id,
                    geom_predicted=dict_evaluated[theme_id],
                    base_formula_field=base_formula_field,
                )
                props[EVALUATION_FIELD_NAME] = Evaluation.TO_CHECK_NO_PREDICTION
                props[PREDICTION_SCORE] = -1
                props[PREDICTION_COUNT] = 0
                if REMARK_FIELD_NAME in props:
                    remarks = props[REMARK_FIELD_NAME]
                else:
                    remarks = []
                remarks.append(ProcessRemark.NO_PREDICTION_ORIGINAL_RETURNED)
                props[REMARK_FIELD_NAME] = remarks
                process_results_evaluated[theme_id][relevant_distance] = unary_union_result_dict({
                    "result": dict_evaluated[theme_id],
                    "properties": props,
                })

        # PART 2: NOT EVALUATED
        relevant_distance = round(0, RELEVANT_DISTANCE_DECIMALS)
        for theme_id, geom in dict_not_evaluated.items():
            process_results_evaluated[theme_id] = {}
            props = self._evaluate(
                id_theme=theme_id,
                geom_predicted=geom,
                base_formula_field=base_formula_field,
            )
            props[EVALUATION_FIELD_NAME] = Evaluation.NOT_EVALUATED
            props[PREDICTION_SCORE] = -1
            if REMARK_FIELD_NAME in props:
                remarks = props[REMARK_FIELD_NAME]
            else:
                remarks = []
            remarks.append(ProcessRemark.NOT_EVALUATED_ORIGINAL_RETURNED)
            props[REMARK_FIELD_NAME] = remarks
            process_results_evaluated[theme_id][relevant_distance] = unary_union_result_dict({
                "result": geom,
                "properties": props,
            })
        return AlignerResult(process_results_evaluated)

    def compare_to_reference(self, geometry: BaseGeometry, with_geom=False):
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
            "reference_source": self.reference_data.source,
            "full": True,
            "area": round(geometry.area, 2),
            "reference_features": {},
            "reference_od": None,
        }

        full_total = True
        last_version_date = None

        ref_intersections = self.reference_data.items.take(
            self.reference_data.tree.query(geometry)
        ).tolist()
        intersected = []
        for key_ref in ref_intersections:
            geom = None
            version_date = None
            geom_reference = self.reference_data[key_ref].geometry
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

    def get_difference_metrics_for_dict_thematic(
        self,
        dict_processresults=None,
        dict_thematic=None,
        diff_metric=DiffMetric.SYMMETRICAL_AREA_CHANGE,
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
            diffs[key] = get_geometry_difference_metrics_from_processresults(
                dict_processresult=dict_processresults[key],
                geom_thematic=dict_thematic[key],
                reference_union=self.reference_data.union,
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
        return self.thematic_data.union

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
        actual_formula = self.compare_to_reference(geom_predicted)
        if actual_formula is None or geom_predicted is None or geom_predicted.is_empty:
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
                area = base_formula["reference_features"][key]["area"]
                if area > 0:
                    diff_percentage_reference_feature = (
                        abs(
                            base_formula["reference_features"][key]["area"]
                            - actual_formula["reference_features"][key]["area"]
                        )
                        * 100
                        / base_formula["reference_features"][key]["area"]
                    )
                else:
                    diff_percentage_reference_feature = 0
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
        # At this moment these are all the check to get a positive EVALUATION. We have to see in future if we add some extra positive EVALUATIONS.
        # fe.
        # * on not-full parcels (comparing all parcels?)
        # * evaluating the outer ring (#102)
        # * ...
        # elif base_formula["full"] == actual_formula["full"] and od_alike:
        #    properties[EVALUATION_FIELD_NAME] = Evaluation.EQUALITY_NON_FULL
        # elif geom_predicted.area >10000: #evaluate only the outer ring
        # pass
        else:
            properties[EVALUATION_FIELD_NAME] = Evaluation.TO_CHECK_NO_PREDICTION
        return properties
