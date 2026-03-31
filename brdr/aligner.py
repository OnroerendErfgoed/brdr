import logging
import os
import uuid
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from concurrent.futures import wait
from datetime import datetime
from typing import Any
from typing import Dict
from typing import Iterable
from typing import List
from typing import Optional
from typing import TYPE_CHECKING
from typing import Union

import numpy as np
from shapely import make_valid, GeometryCollection
from shapely import to_geojson
from shapely.geometry.base import BaseGeometry
from shapely.prepared import prep

from brdr import __version__
from brdr.configs import AlignerConfig
from brdr.configs import ProcessorConfig
from brdr.constants import (
    AREA_CHANGE,
    METADATA_FIELD_NAME,
    MAX_REFERENCE_BUFFER,
    OBSERVATION_FIELD_NAME,
)
from brdr.constants import AREA_PERCENTAGE_CHANGE
from brdr.constants import DATE_FORMAT
from brdr.constants import DEFAULT_CRS
from brdr.constants import EVALUATION_FIELD_NAME
from brdr.constants import ID_THEME_FIELD_NAME
from brdr.constants import LAST_VERSION_DATE
from brdr.constants import LENGTH_CHANGE
from brdr.constants import LENGTH_PERCENTAGE_CHANGE
from brdr.constants import NR_CALCULATION_FIELD_NAME
from brdr.constants import PREDICTION_SCORE
from brdr.constants import RELEVANT_DISTANCE_FIELD_NAME
from brdr.constants import REMARK_FIELD_NAME
from brdr.constants import SYMMETRICAL_AREA_CHANGE
from brdr.constants import SYMMETRICAL_AREA_PERCENTAGE_CHANGE
from brdr.constants import VERSION_DATE
from brdr.enums import AlignerResultType
from brdr.enums import DiffMetric
from brdr.enums import Evaluation
from brdr.enums import FullReferenceStrategy
from brdr.enums import ProcessRemark
from brdr.evaluator import AlignerEvaluator
from brdr.evaluator import BaseEvaluator
from brdr.feature_data import AlignerFeatureCollection
from brdr.geometry_utils import buffer_neg
from brdr.geometry_utils import buffer_pos
from brdr.geometry_utils import safe_difference
from brdr.geometry_utils import safe_intersection
from brdr.geometry_utils import safe_unary_union
from brdr.geometry_utils import to_crs
from brdr.loader import Loader
from brdr.logger import Logger
from brdr.metadata import (
    get_metadata_observations_from_process_result as _get_metadata_observations_core,
)
from brdr.metadata import (
    reverse_metadata_observations_to_brdr_observation as _reverse_metadata_observations_core,
)
from brdr.predictor import AlignerPredictor
from brdr.predictor import BasePredictor
from brdr.processor import AlignerGeometryProcessor
from brdr.processor import BaseProcessor
from brdr.typings import InputId
from brdr.typings import ProcessResult
from brdr.utils import get_geodataframe_from_process_results
from brdr.utils import get_geojsons_from_process_results
from brdr.utils import get_geometry_difference_metrics_from_processresult
from brdr.utils import get_geometry_difference_metrics_from_processresults
from brdr.utils import urn_from_geom
from brdr.utils import write_featurecollection_to_geojson

###################

if TYPE_CHECKING:
    from brdr.aligner import Aligner


class AlignerResult:
    """
    Stores and processes the results from the alignment operation.

    This class manages the lifecycle of alignment results, including enrichment with
    geometric metrics (like area or length change), filtering by result type,
    and exporting to GeoJSON format.

    Attributes
    ----------
    results : Dict[InputId, Dict[float, Optional[ProcessResult]]]
        A nested dictionary where the outer key is the theme ID and the inner key
        is the relevant distance used during processing.
    metadata : Optional[Any]
        Optional metadata about the overall alignment operation.

    Notes
    -----
    The class acts as a post-processor. It takes the raw output from the geometric
    engines and calculates semantic differences before exporting.

    ```{mermaid}
    graph TD
        Raw[Raw Process Results] --> AR[AlignerResult]
        AR --> |get_results| Enriched[Enriched with Area/Length Metrics]
        Enriched --> |get_results_as_geojson| GeoJSON[GeoJSON FeatureCollection]
        GeoJSON --> |save_results| Disk[Local .geojson files]
    ```

    Examples
    --------
    >>> # Initializing and saving results
    >>> result_obj = AlignerResult(raw_results)
    >>> result_obj.save_results(aligner, path="./output_folder")
    """

    metadata: Optional[Any] = None

    def __init__(
        self,
        process_results: Dict[InputId, Dict[float, Optional[ProcessResult]]],
    ):
        """
        Initializes an AlignerResult instance.

        Parameters
        ----------
        process_results : Dict[InputId, Dict[float, Optional[ProcessResult]]]
            A nested dictionary containing raw results.
            Structure: `{theme_id: {distance: result_object}}`.
        """
        self.results = process_results

    def get_results(
        self,
        aligner: "Aligner",
        result_type: AlignerResultType = AlignerResultType.PROCESSRESULTS,
    ) -> Dict[InputId, Dict[float, Optional[ProcessResult]]]:
        """
        Retrieves results enriched with geometric metrics and filtered by type.

        This method mutates the stored `ProcessResult` objects by adding calculated
        difference metrics (e.g., symmetrical area change) to their properties.

        Parameters
        ----------
        aligner : Aligner
            The 'Aligner' object used to access thematic data and comparison logic.
        result_type : AlignerResultType
            Filters the results based on status. Defaults to PROCESSRESULTS.

        Returns
        -------
        Dict[InputId, Dict[float, Optional[ProcessResult]]]
            A dictionary of filtered and enriched results.

        Raises
        ------
        ValueError
            If an unsupported `result_type` is provided.
        """
        for theme_id, results_dict in self.results.items():
            original_feature = aligner.thematic_data.features.get(theme_id)
            if original_feature is None:
                continue

            original_geometry = original_feature.geometry
            nr_calculations = len(results_dict)

            try:
                original_geometry_length = len(original_geometry.geoms)
            except (AttributeError, TypeError):
                original_geometry_length = 1

            for relevant_distance, process_result in results_dict.items():
                if process_result is None:
                    continue

                properties = process_result["properties"]
                properties[ID_THEME_FIELD_NAME] = theme_id
                properties[NR_CALCULATION_FIELD_NAME] = nr_calculations
                properties[RELEVANT_DISTANCE_FIELD_NAME] = relevant_distance

                # Add Geometric Difference Metrics
                metrics_to_calc = [
                    (SYMMETRICAL_AREA_CHANGE, DiffMetric.SYMMETRICAL_AREA_CHANGE),
                    (
                        SYMMETRICAL_AREA_PERCENTAGE_CHANGE,
                        DiffMetric.SYMMETRICAL_AREA_PERCENTAGE_CHANGE,
                    ),
                    (AREA_CHANGE, DiffMetric.AREA_CHANGE),
                    (AREA_PERCENTAGE_CHANGE, DiffMetric.AREA_PERCENTAGE_CHANGE),
                    (LENGTH_CHANGE, DiffMetric.LENGTH_CHANGE),
                    (LENGTH_PERCENTAGE_CHANGE, DiffMetric.LENGTH_PERCENTAGE_CHANGE),
                ]

                for prop_key, metric_enum in metrics_to_calc:
                    properties[prop_key] = (
                        get_geometry_difference_metrics_from_processresult(
                            process_result, original_geometry, None, metric_enum
                        )
                    )

                # Check for Geometry Count Change
                resulting_geom = process_result["result"]
                try:
                    resulting_geometry_length = len(resulting_geom.geoms)
                except (AttributeError, TypeError):
                    resulting_geometry_length = 1

                if original_geometry_length != resulting_geometry_length:
                    remark = ProcessRemark.CHANGED_AMOUNT_GEOMETRIES
                    remarks = properties.get(REMARK_FIELD_NAME, [])
                    remarks.append(remark)
                    properties[REMARK_FIELD_NAME] = remarks

        if result_type == AlignerResultType.PROCESSRESULTS:
            return self.results
        elif result_type == AlignerResultType.PREDICTIONS:
            return {
                tid: {
                    rd: res
                    for rd, res in rd_dict.items()
                    if res and PREDICTION_SCORE in res["properties"]
                }
                for tid, rd_dict in self.results.items()
            }
        elif result_type == AlignerResultType.EVALUATED_PREDICTIONS:
            return {
                tid: {
                    rd: res
                    for rd, res in rd_dict.items()
                    if res
                    and EVALUATION_FIELD_NAME
                    in res[
                        "properties"
                    ]  # and res["properties"][EVALUATION_FIELD_NAME]!=Evaluation.NOT_EVALUATED
                }
                for tid, rd_dict in self.results.items()
            }
        else:
            raise ValueError(f"Unknown result type: {result_type}")

    def get_results_as_geojson(
        self,
        aligner: "Aligner",
        result_type: AlignerResultType = AlignerResultType.PROCESSRESULTS,
        add_metadata: bool = False,
        add_original_attributes: bool = False,
    ) -> Dict[str, Any]:
        """
        Converts the results into a GeoJSON FeatureCollection format.

        Parameters
        ----------
        aligner : Aligner
            The 'Aligner' object providing CRS and thematic metadata.
        result_type : AlignerResultType
            The type of results to export. Defaults to PROCESSRESULTS.
        add_metadata : bool
            If True, includes the descriptive comparison observation in properties.
        add_original_attributes : bool
            If True, includes original thematic attributes in the output.

        Returns
        -------
        Dict[str, Any]
            A dictionary representing a GeoJSON FeatureCollection.

        Raises
        ------
        ValueError
            If `self.results` is empty or None.
        """
        if not self.results:
            raise ValueError("Empty results: No calculated results to export.")

        results = self.get_results(aligner=aligner, result_type=result_type)
        prop_dictionary = defaultdict(dict)

        for theme_id, results_dict in results.items():
            for relevant_distance, process_result in results_dict.items():
                if process_result is None:
                    continue

                prop_dictionary[theme_id][relevant_distance] = {}
                feature = aligner.thematic_data.features.get(theme_id)

                if add_original_attributes and feature and feature.properties:
                    prop_dictionary[theme_id][relevant_distance].update(
                        feature.properties
                    )

                if add_metadata:
                    metadata_result = process_result.get("metadata", None)
                    if metadata_result is None:
                        logging.debug("metadata not available")
                        continue
                    prop_dictionary[theme_id][relevant_distance][
                        METADATA_FIELD_NAME
                    ] = metadata_result

        return get_geojsons_from_process_results(
            results,
            crs=aligner.crs,
            id_field=aligner.thematic_data.id_fieldname,
            series_prop_dict=prop_dictionary,
        )

    def get_results_as_geodataframe(
        self,
        aligner: "Aligner",
        result_type: AlignerResultType = AlignerResultType.PROCESSRESULTS,
        add_metadata: bool = False,
        add_original_attributes: bool = False,
    ):
        """
        Converts the results into a flat GeoDataFrame.

        The resulting GeoDataFrame contains one row per
        (theme_id, relevant_distance), geometry columns (result, result_diff, ...),
        and flattened dictionary columns (e.g. properties__*).
        """
        if not self.results:
            raise ValueError("Empty results: No calculated results to export.")

        results = self.get_results(aligner=aligner, result_type=result_type)
        prop_dictionary = defaultdict(dict)

        for theme_id, results_dict in results.items():
            for relevant_distance, process_result in results_dict.items():
                if process_result is None:
                    continue

                prop_dictionary[theme_id][relevant_distance] = {}
                feature = aligner.thematic_data.features.get(theme_id)

                if add_original_attributes and feature and feature.properties:
                    prop_dictionary[theme_id][relevant_distance].update(
                        feature.properties
                    )

                if add_metadata:
                    metadata_result = process_result.get("metadata", None)
                    if metadata_result is None:
                        logging.debug("metadata not available")
                        continue
                    prop_dictionary[theme_id][relevant_distance][
                        METADATA_FIELD_NAME
                    ] = metadata_result

        return get_geodataframe_from_process_results(
            results,
            crs=aligner.crs,
            id_field=aligner.thematic_data.id_fieldname,
            series_prop_dict=prop_dictionary,
            preferred_geometry_column="result",
        )

    def save_results(
        self,
        aligner: "Aligner",
        path: str,
        result_type: AlignerResultType = AlignerResultType.PROCESSRESULTS,
        add_metadata: bool = False,
        add_original_attributes: bool = False,
    ) -> None:
        """
        Exports the results as multiple GeoJSON files to a directory.

        Creates separate files for original results, differences, and specific
        area changes (added/removed).

        Parameters
        ----------
        aligner : Aligner
            The 'Aligner' object for spatial context.
        path : str
            Target directory path where files will be created.
        result_type : AlignerResultType
            Type of results to export. Defaults to PROCESSRESULTS.
        add_metadata : bool
            Whether to include alignment observations in the output.
        add_original_attributes : bool
            Whether to include original feature attributes.

        Notes
        -----
        The output files follow the naming convention: `{result_type}_{name}.geojson`.
        """
        fcs = self.get_results_as_geojson(
            aligner=aligner,
            result_type=result_type,
            add_metadata=add_metadata,
            add_original_attributes=add_original_attributes,
        )
        for name, fc in fcs.items():
            file_name = f"{result_type.value}_{name}.geojson"
            write_featurecollection_to_geojson(os.path.join(path, file_name), fc)


def _get_metadata_observations_from_process_result(
    processResult: ProcessResult,
    reference_lookup: Dict[any, any],
) -> List[Dict]:
    return _get_metadata_observations_core(processResult, reference_lookup)


def _reverse_metadata_observations_to_brdr_observation(metadata: Dict) -> Dict:
    """
    Reconstruct a BRDR observation dictionary from a list of Linked Data observations.

    This function performs a reverse transformation by aggregating individual
    SOSA-based observations back into a nested BRDR-structured dictionary.
    It maps RDF-style 'observed_property' URIs back to their original JSON keys
    and reconstructs the reference features and actuation metadata.

    Parameters
    ----------
    metadata : List[Dict]
        A list containing observation dictionaries. Expected to have a key
        "observations" which holds a list of dictionaries following the
        SOSA (Sensor, Observation, Sample, and Actuator) ontology.

    Returns
    -------
    Dict
        A dictionary representing the reconstructed BRDR observation,
        including 'reference_features' and global observation metrics.

    Raises
    ------
    ValueError
        If the reconstructed dictionary does not pass the `is_brdr_observation`
        validation check.

    Notes
    -----
    The function handles the de-duplication of reference geometries by tracking
    unique IDs in the 'used' field of the observations. It distinguishes
    between root-level properties and reference-specific properties by
    checking for the presence of a 'ref_id' that differs from the primary
    result interest ID.
    """
    return _reverse_metadata_observations_core(metadata)


def aligner_metadata_decorator(f):
    def inner_func(thematic_id, geometry, relevant_distance, aligner, *args, **kwargs):
        process_result: ProcessResult = f(
            thematic_id, geometry, relevant_distance, aligner, *args, **kwargs
        )
        if aligner.add_observations:
            process_result["observations"] = aligner.compare_to_reference(
                process_result.get("result")
            )
            # add observation properties to the properties
            observation_props = aligner.get_observation_properties(process_result)
            props = process_result["properties"]
            props[OBSERVATION_FIELD_NAME] = process_result[
                "observations"
            ]  # adding the raw brdr_observation?!
            props.update(observation_props)
            process_result["properties"] = props

        if aligner.log_metadata:
            # generate uuid for actuation
            actuation_id = uuid.uuid4()
            processor_id = aligner.processor.processor_id.value
            processor_name = type(aligner.processor).__name__
            reference_data = aligner.reference_data
            reference_intersections_ids = reference_data.items.take(
                reference_data.tree.query(buffer_pos(geometry, MAX_REFERENCE_BUFFER))
            ).tolist()  # TODO possible to optimize?
            reference_geometries = []
            for ref_id in reference_intersections_ids:
                feature = reference_data.features[ref_id]
                feat_dict = {
                    # "id": feature.data_id,
                    "id": feature.brdr_id,
                    "type": f"geo:{feature.geometry.geom_type}",
                    "version_date": reference_data.source.get(VERSION_DATE, ""),
                    "derived_from": {
                        "id": feature.data_id,
                        "type": "geo:Feature",
                        "source": reference_data.source.get("source_url", ""),
                    },
                }
                reference_geometries.append(feat_dict)

            thematic_feature = aligner.thematic_data.features[thematic_id]
            feature_of_interest_id = thematic_feature.brdr_id
            result_urn = urn_from_geom(process_result["result"])
            actuation_time = datetime.now().strftime("%Y-%m-%dT%H:%M:%S%z")
            process_result["metadata"] = {}
            process_result["metadata"]["actuation"] = {
                "id": actuation_id.urn,
                "type": "sosa:Actuation",
                "reference_geometries": reference_geometries,
                "changes": "geo:hasGeometry",
                "sosa:hasFeatureOfInterest": {"id": feature_of_interest_id},
                "result": result_urn,
                "result_time": actuation_time,  # TODO _CHECK
                "procedure": {
                    "id": processor_id,
                    "implementedBy": processor_name,
                    "type": "sosa:Procedure",
                    "ssn:hasInput": [
                        {
                            "id": "brdr:relevant_distance",
                            "type": "ssn:Input",
                            "input_value": {
                                "type": "xsd:integer",
                                "value": relevant_distance,
                            },
                        },
                    ],
                },
            }
            if process_result["observations"]:
                ref_lookup = reference_data.reference_lookup
                process_result["metadata"]["observations"] = (
                    _get_metadata_observations_from_process_result(
                        process_result, ref_lookup
                    )
                )

        return process_result

    return inner_func


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
    add_observations: bool
        If True, process result observations will be computed by default.
    processor : BaseProcessor or AlignerGeometryProcessor
        The geometric processor used for alignment calculations.
    predictor : BasePredictor or AlignerPredictor
        The predictor strategy used to assign prediction scores on process results.
    evaluator : BaseEvaluator or AlignerEvaluator
        The evaluator strategy used to evaluate predicted candidates.
    correction_distance : float
        Distance used in buffer operations to remove slivers (technical correction).
    mitre_limit : int
        Parameter for the buffer operation to control the maximum length of join corners.
    max_workers : int, optional
        The maximum number of workers for parallel execution (ThreadPoolExecutor).
    crs : str
        The Coordinate Reference System (CRS) being used (bv. 'EPSG:31370').
    name_thematic_id : str
        Name of the identifier field for thematic data.
    diff_metric : DiffMetric
        The metric used to measure differences between geometries (bv. area change).
    reference_data : AlignerFeatureCollection, optional
        Loaded collection of reference features.
    thematic_data : AlignerFeatureCollection, optional
        Loaded collection of thematic features.

    Notes
    -----
    The Aligner acts as the central orchestrator of the `brdr` package. It connects
    data loaders with geometric processors.

    ```{mermaid}
    graph LR
        T[Thematic Data] --> Aligner
        R[Reference Data] --> Aligner
        Aligner --> P{Processor}
        P --> Res[AlignerResult]
    ```

    Examples
    --------
    >>> from brdr.aligner import Aligner
    >>> aligner = Aligner(crs="EPSG:31370")
    >>> aligner.load_thematic_data(...)
    >>> aligner.load_reference_data(...)
    >>> results = aligner.process(relevant_distances=[0, 0.5, 1.0])
    """

    def __init__(
        self,
        *,
        processor: Optional[BaseProcessor] = None,
        predictor: Optional[BasePredictor] = None,
        evaluator: Optional[BaseEvaluator] = None,
        crs: str = DEFAULT_CRS,
        config: Optional[AlignerConfig] = None,
        feedback: Any = None,
    ):
        """
        Initializes the Aligner object.

        Parameters
        ----------
        processor : BaseProcessor, optional
            The geometric processor instance. If None, AlignerGeometryProcessor is used.
        predictor : BasePredictor, optional
            The prediction strategy instance. If None, AlignerPredictor is used.
        evaluator : BaseEvaluator, optional
            The evaluation strategy instance. If None, AlignerEvaluator is used.
        crs : str, optional
            Coordinate Reference System (CRS) of the data.
            Expected to be a projected CRS with units in meters.
            Defaults to DEFAULT_CRS (EPSG:31370).
        config : AlignerConfig, optional
            Configuration object containing parameters like correction_distance,
            diff_metric, and max_workers.
        feedback : object, optional
            Feedback object used for logging (e.g., in a QGIS environment).
            Defaults to None.

        Notes
        -----
        If a `config` object is provided, its values will override the default
        settings for internal attributes like `correction_distance` and `mitre_limit`.
        """
        if config is None:
            config = AlignerConfig()
        self.config = config
        self.logger = Logger(feedback)
        self.log_metadata = config.log_metadata
        self.add_observations = config.add_observations
        self.processor = (
            processor
            if processor
            else AlignerGeometryProcessor(ProcessorConfig(), feedback)
        )
        self.predictor = predictor if predictor else AlignerPredictor(feedback)
        self.evaluator = evaluator if evaluator else AlignerEvaluator(feedback)
        self.correction_distance = config.correction_distance
        self.mitre_limit = config.mitre_limit
        self.max_workers = config.max_workers

        # PROCESSING DEFAULTS (Internal state variables)

        # The CRS is the working CRS for all calculations (assumed to be projected/in meters)
        self.crs = to_crs(crs)

        self.diff_metric = config.diff_metric
        self.logger.feedback_info("Aligner initialized")

        self.reference_data: AlignerFeatureCollection | None = None
        self.thematic_data: AlignerFeatureCollection | None = None

    ##########LOADERS##########################
    ###########################################

    def load_thematic_data(self, loader: Loader):
        """
        Loads the thematic features into the aligner using a specific loader.

        This method executes the loader's data retrieval logic and ensures that
        the resulting feature collection is tagged with the Aligner's CRS.

        Parameters
        ----------
        loader : Loader
            An instance of a Loader class (e.g., GeoJsonLoader or WFSReferenceLoader)
            that implements the `load_data` interface.

        Notes
        -----
        The method automatically synchronizes the CRS of the loaded data with
        the `self.crs` attribute of the Aligner instance.

        Examples
        --------
        >>> from brdr.loader import GeoJsonLoader
        >>> loader = GeoJsonLoader(path="data/themes.json")
        >>> aligner.load_thematic_data(loader)
        """

        self.thematic_data = loader.load_data()
        self.thematic_data.crs = self.crs

    def load_reference_data(self, loader: Loader):
        """
        Loads the reference features into the aligner and prepares them for processing.

        This method retrieves data via the provided loader, synchronizes the CRS,
        and marks the feature collection as a reference dataset to enable
        spatial indexing and comparison logic.

        Parameters
        ----------
        loader : Loader
            An instance of a Loader class (e.g., WFSReferenceLoader) that implements
            the `load_data` interface.

        Notes
        -----
        Setting `is_reference = True` on the resulting dataset is essential for
        the internal alignment logic, as it distinguishes the 'anchor' geometries
        from the thematic geometries that need to be shifted.

        Examples
        --------
        >>> from brdr.loader import WFSReferenceLoader
        >>> ref_loader = WFSReferenceLoader(url="https://geoserver.com/wfs")
        >>> aligner.load_reference_data(ref_loader)
        >>> print(aligner.reference_data.is_reference)
        True
        """

        self.reference_data = loader.load_data()
        self.reference_data.crs = self.crs
        self.reference_data.is_reference = True

    def process(
        self,
        relevant_distances: Iterable[float] = None,
        *,
        thematic_ids: List[InputId] = None,
        max_workers: int = None,
    ) -> AlignerResult:
        """
        Executes the alignment process across multiple relevant distances.

        This method iterates through the specified thematic features and calculates
        the alignment for each provided 'relevant distance'. It can utilize
        parallel processing via a ThreadPoolExecutor to speed up calculations.

        Parameters
        ----------
        relevant_distances : Iterable[float]
            A series of distances (in meters) to use for the alignment logic.
            This parameter is mandatory.
        thematic_ids : List[InputId], optional
            A specific list of IDs to process. If None, all thematic features
            currently loaded in the aligner will be processed.
        max_workers : int, optional
            The number of threads for parallel execution. If -1, execution is
            serial. If None, the Aligner's default `max_workers` is used.

        Returns
        -------
        [AlignerResult][]
            An object containing the structured results, accessible by thematic ID
            and relevant distance.

        Raises
        ------
        ValueError
            If `relevant_distances` is None or empty.
            If any provided `thematic_ids` are not found in the loaded thematic data.

        Notes
        -----
        The processing flow involves distributing geometry-distance pairs across
        available worker threads:

        ```{mermaid}
        graph TD
            Start[Process Call] --> Check{Valid IDs?}
            Check -- Yes --> Parallel{max_workers > 0?}
            Parallel -- Yes --> Workers[ThreadPoolExecutor]
            Parallel -- No --> Serial[Single Thread]
            Workers --> Proc[Geometry Processor]
            Serial --> Proc
            Proc --> Result[AlignerResult]
        ```

        Examples
        --------
        >>> aligner.process(relevant_distances=[0, 0.5, 1.0], max_workers=4)
        <brdr.aligner.AlignerResult object at ...>
        """
        if relevant_distances is None:
            raise ValueError("provide at least 1 relevant distance")
        if thematic_ids is None:
            thematic_ids = self.thematic_data.features.keys()
        if any(
            id_to_process not in self.thematic_data.features.keys()
            for id_to_process in thematic_ids
        ):
            raise ValueError("not all ids are found in the thematic data")
        if max_workers is None:
            max_workers = self.max_workers

        self.logger.feedback_debug("Process series" + str(relevant_distances))

        process_results: dict[InputId, dict[float, ProcessResult | None]] = {}
        futures = {}

        @aligner_metadata_decorator
        def process_geom_for_rd(thematic_id, geometry, relevant_distance, aligner=self):
            return aligner.processor.process(
                correction_distance=aligner.correction_distance,
                reference_data=aligner.reference_data,
                mitre_limit=aligner.mitre_limit,
                input_geometry=geometry,
                relevant_distance=relevant_distance,
                thematic_data=aligner.thematic_data,
            )

        def run_process(executor: ThreadPoolExecutor = None):
            for thematic_id in thematic_ids:
                geom = self.thematic_data.features.get(thematic_id).geometry
                self.logger.feedback_debug(
                    f"thematic id {str(thematic_id)} processed with "
                    f"relevant distances (m) [{str(relevant_distances)}]"
                )
                process_results[thematic_id] = {}
                for rd in relevant_distances:
                    try:
                        fn = process_geom_for_rd
                        if executor:
                            futures[(thematic_id, rd)] = executor.submit(
                                fn, thematic_id, geom, rd, self
                            )
                        else:
                            process_results[thematic_id][rd] = fn(
                                thematic_id, geom, rd, self
                            )
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
        relevant_distances: Optional[Union[List[float], np.ndarray]] = None,
        *,
        thematic_ids: Optional[List[InputId]] = None,
        diff_metric: Optional[DiffMetric] = None,
    ) -> AlignerResult:
        """
        Predicts the 'most interesting' relevant distances for changes in thematic
        elements.

        This method analyzes the stability of geometric differences across a range
        of distances. It identifies "breakpoints" where changes occur and
        "zero-streaks" where the geometry remains stable. Based on this, a
        prediction score is assigned to specific results.

        Parameters
        ----------
        relevant_distances : List[float] or np.ndarray, optional
            A series of distances (in meters) to be analyzed.
            Defaults to a range from 0.0 to 3.0 meters with steps of 0.1m.
        thematic_ids : List[InputId], optional
            Specific thematic IDs to process. If None, all loaded thematic
            geometries are used.
        diff_metric : DiffMetric, optional
            The metric used to determine differences (e.g., area change).
            If None, the Aligner's default `diff_metric` is used.

        Returns
        -------
        AlignerResult
            A result object containing the process results enriched with
            stability properties and prediction scores.

        Raises
        ------
        ValueError
            If provided `thematic_ids` are not found in the loaded data.

        Notes
        -----
        The prediction logic follows a four-step process:
        1. **Series Processing**: Calculates results for all specified distances.
        2. **Metric Calculation**: Computes the difference between original and
           result for each step.
        3. **Stability Analysis**: Identifies streaks where the change is minimal
           (zero-streaks) or where the metric flips sign (breakpoints).
        4. **Scoring**: Assigns a `PREDICTION_SCORE` to results that fall within
           a significant stability window.

        ```{mermaid}
        graph TD
            A[Input Distances] --> B[Process Series]
            B --> C[Calculate Diffs]
            C --> D[Stability Analysis]
            D --> E{Stable?}
            E -- Yes --> F[Assign Prediction Score]
            E -- No --> G[Mark unstable]
            F --> H[AlignerResult]
            G --> H
        ```

        Examples
        --------
        >>> aligner = Aligner(crs="EPSG:31370")
        >>> # Predict using default distance range
        >>> prediction_results = aligner.predict(thematic_ids=["id_1", "id_2"])
        """
        return self.predictor.predict(
            aligner=self,
            relevant_distances=relevant_distances,
            thematic_ids=thematic_ids,
            diff_metric=diff_metric,
        )

    def evaluate(
        self,
        relevant_distances: Optional[Iterable[float]] = None,
        *,
        thematic_ids: Optional[List[InputId]] = None,
        metadata_field: str = METADATA_FIELD_NAME,
        full_reference_strategy: FullReferenceStrategy = FullReferenceStrategy.NO_FULL_REFERENCE,
        max_predictions: int = -1,
        multi_to_best_prediction: bool = True,
    ) -> AlignerResult:
        """
        Evaluates input geometries against predictions using observation matching
        and selection strategies.

        This method identifies the best candidate geometries for alignment by
        comparing predicted geometries against the original attributes and
        reference data. It tags results with evaluation labels (e.g.,
        PREDICTION_UNIQUE) to facilitate decision-making.

        Parameters
        ----------
        relevant_distances : Iterable[float], optional
            Distances to evaluate. Defaults to 0.0m to 3.0m (step 0.1m).
        thematic_ids : List[InputId], optional
            List of IDs to evaluate. If None, all loaded thematic features
            are processed. Features not in this list are marked as NOT_EVALUATED.
        metadata_field : str, optional
            The field name containing the metadata of the input geometry.
            Defaults to METADATA_FIELD_NAME.
        full_reference_strategy : FullReferenceStrategy, optional
            Strategy to prioritize predictions based on their 'fullness'
            relative to reference data. Defaults to NO_FULL_REFERENCE.
        max_predictions : int, optional
            Maximum number of predictions to return per feature.
            -1 returns all candidates. Defaults to -1.
        multi_to_best_prediction : bool, optional
            If True and `max_predictions=1`, returns the candidate with the
            highest score. If False, returns the original geometry when
            multiple candidates exist. Defaults to True.

        Returns
        -------
        AlignerResult
            An object containing evaluated results with detailed metadata and
            evaluation status tags.

        Raises
        ------
        ValueError
            If provided `thematic_ids` are not present in the thematic data.

        Notes
        -----
        The evaluation logic follows a specific priority tree:
        1. **observation Match**: If a prediction matches the base observation exactly,
           it is prioritized (score 100).
        2. **Fullness**: If a strategy is set, 'full' predictions get a score boost.
        3. **Scoring**: Candidates are ranked by their prediction score.
        4. **Selection**: Based on `max_predictions`, the best candidates are selected.



        ```{mermaid}
        graph TD
            Start[Start Evaluation] --> Pred[Generate Predictions]
            Pred --> Match{observation Match?}
            Match -- Yes --> HighScore[Prioritize & Score 100]
            Match -- No --> Full{Full Reference?}
            Full -- Yes --> Strategy[Apply Fullness Strategy]
            Full -- No --> Rank[Rank by Prediction Score]
            Strategy --> Rank
            HighScore --> Rank
            Rank --> Select{Max Predictions?}
            Select --> Final[Return AlignerResult]
        ```

        Examples
        --------
        >>> aligner = Aligner()
        >>> # Evaluate with a limit of 1 best prediction per feature
        >>> results = aligner.evaluate(thematic_ids=["id_01"],full_reference_strategy=FullReferenceStrategy.ONLY_FULL_REFERENCE,max_predictions=1)
        """
        return self.evaluator.evaluate(
            aligner=self,
            relevant_distances=relevant_distances,
            thematic_ids=thematic_ids,
            metadata_field=metadata_field,
            full_reference_strategy=full_reference_strategy,
            max_predictions=max_predictions,
            multi_to_best_prediction=multi_to_best_prediction,
        )

    def _update_evaluation_with_original(
        self,
        metadata_field: str,
        original_geometry: BaseGeometry,
        process_results_evaluated,
        theme_id: str | int,
        evaluation: Evaluation,
    ) -> Any:
        return self.evaluator.update_evaluation_with_original(
            aligner=self,
            metadata_field=metadata_field,
            original_geometry=original_geometry,
            process_results_evaluated=process_results_evaluated,
            theme_id=theme_id,
            evaluation=evaluation,
        )

    def compare_to_reference(
        self, geometry: BaseGeometry, with_geom: bool = False
    ) -> Dict[str, Any]:
        """
        Calculates observation-related information based on the input geometry.

        This method performs a spatial analysis to determine how much of the input
        geometry is covered by reference features and identifies which specific
        reference elements are involved. It also detects "OD" (Onbekend Terrein/Open
        Space) areas that are not covered by any reference data.

        Parameters
        ----------
        geometry : BaseGeometry
            The input geometry to be analyzed against the reference dataset.
        with_geom : bool, optional
            If True, includes the GeoJSON representation of each intersection
            and the OD-area in the output. Defaults to False.

        Returns
        -------
        Dict[str, Any]
            A dictionary containing the alignment analysis:

            - **alignment_date** (str): Timestamp of the calculation.
            - **brdr_version** (str): Version of the package used.
            - **full** (bool): True if the geometry is entirely composed of
              complete reference features.
            - **area** (float): Total area of the input geometry.
            - **last_version_date** (str, optional): The most recent version date
              found among the intersected reference features.
            - **reference_features** (dict): A mapping of reference IDs to:
                - *full* (bool): True if the reference feature is fully contained.
                - *area* (float): Area of the intersection.
                - *percentage* (float): Coverage percentage relative to the reference.
                - *geometry* (str, optional): GeoJSON string (if with_geom is True).
            - **reference_od** (dict, optional): Description of areas not covered
              by reference features, containing 'area' and optionally 'geometry'.

        Notes
        -----
        The method uses a spatial index (R-tree) to efficiently find intersecting
        reference features. A correction distance is applied to the "OD" calculation
        to filter out insignificant slivers.



        ```{mermaid}
        graph TD
            In[Input Geometry] --> Query[Spatial Index Query]
            Query --> Intersect[Calculate Intersections]
            Intersect --> Stats[Compute Area & %]
            In --> Diff[Difference with Union of Refs]
            Diff --> OD[Identify Unknown Area / OD]
            Stats --> Dict[Final observation Dictionary]
            OD --> Dict
        ```

        Examples
        --------
        >>> info = aligner.compare_to_reference(my_geom, with_geom=True)
        >>> print(info["full"])
        >>> print(info["reference_features"].keys())  # IDs of intersected refs
        """
        def _count_points(geom: BaseGeometry) -> float:
            if geom is None or geom.is_empty:
                return 0.0
            if geom.geom_type == "Point":
                return 1.0
            if geom.geom_type == "MultiPoint":
                return float(len(geom.geoms))
            if geom.geom_type == "GeometryCollection":
                return float(sum(_count_points(g) for g in geom.geoms))
            return 0.0

        def _measure_type_for_geom(geom: BaseGeometry) -> str:
            if geom.geom_type in {"Point", "MultiPoint"}:
                return "count"
            if geom.geom_type in {"LineString", "MultiLineString"}:
                return "length"
            return "area"

        def _reference_metric_for_pair(
            thematic_geom: BaseGeometry,
            reference_geom: BaseGeometry,
            intersection_geom: BaseGeometry,
        ) -> tuple[str, float, float]:
            thematic_type = _measure_type_for_geom(thematic_geom)
            reference_type = _measure_type_for_geom(reference_geom)

            # Polygon thematic geometry:
            # - polygon refs: area overlap over reference area (existing behavior)
            # - line refs: boundary coverage length over reference length
            # - point refs: intersecting points over reference point count
            if thematic_type == "area":
                if reference_type == "area":
                    return (
                        "area",
                        float(intersection_geom.area),
                        float(reference_geom.area),
                    )
                if reference_type == "length":
                    return (
                        "length",
                        float(intersection_geom.length),
                        float(reference_geom.length),
                    )
                return (
                    "count",
                    _count_points(intersection_geom),
                    _count_points(reference_geom),
                )

            # Line thematic geometry:
            # - line/area refs: overlapping line length over thematic length
            # - point refs: touched reference points over reference point count
            if thematic_type == "length":
                if reference_type == "count":
                    return (
                        "count",
                        _count_points(intersection_geom),
                        _count_points(reference_geom),
                    )
                return (
                    "length",
                    float(intersection_geom.length),
                    float(thematic_geom.length),
                )

            # Point thematic geometry:
            # - always point counts over thematic point count
            return (
                "count",
                _count_points(intersection_geom),
                _count_points(thematic_geom),
            )

        def _measure_value(geom: BaseGeometry, measure_type: str) -> float:
            if geom is None or geom.is_empty:
                return 0.0
            if measure_type == "count":
                return _count_points(geom)
            if measure_type == "length":
                return float(geom.length)
            return float(geom.area)

        measure_type = _measure_type_for_geom(geometry)
        thematic_measure = _measure_value(geometry, measure_type)

        dict_observation = {
            "alignment_date": datetime.now().strftime(DATE_FORMAT),
            "brdr_version": str(__version__),
            "reference_source": self.reference_data.source,
            "full": True,
            "area": round(geometry.area, 2),
            "length": round(geometry.length, 2),
            "count": round(_count_points(geometry), 2),
            "measure_type": measure_type,
            "reference_features": {},
            "reference_od": None,
        }

        full_total = True
        last_version_date = None
        prepared_geometry = prep(geometry)

        ref_intersections = self.reference_data.items.take(
            self.reference_data.tree.query(geometry)
        ).tolist()
        intersected = []
        for key_ref in ref_intersections:
            geom = None
            version_date = None
            geom_reference = self.reference_data[key_ref].geometry
            if not prepared_geometry.intersects(geom_reference):
                continue
            geom_intersection = make_valid(safe_intersection(geometry, geom_reference))
            if geom_intersection.is_empty or geom_intersection is None:
                continue
            intersected.append(geom_intersection)

            (
                feature_measure_type,
                intersection_measure,
                denominator,
            ) = _reference_metric_for_pair(
                thematic_geom=geometry,
                reference_geom=geom_reference,
                intersection_geom=geom_intersection,
            )
            if denominator > 0:
                perc = round(intersection_measure * 100 / denominator, 2)
            else:
                perc = 0
            if perc < 0.01:
                continue
            # Add a last_version_date if available in properties
            reference_feature = self.reference_data.features.get(key_ref)
            if (
                not reference_feature is None
                and VERSION_DATE in reference_feature.properties
            ):
                str_version_date = reference_feature.properties[VERSION_DATE]
                version_date = datetime.strptime(str_version_date, DATE_FORMAT)
                if last_version_date is None and version_date is not None:
                    last_version_date = version_date
                if version_date is not None and version_date > last_version_date:
                    last_version_date = version_date

            if perc > 99.99:
                full = True
                area = round(
                    geom_reference.area
                    if measure_type == "area"
                    else geom_intersection.area,
                    2,
                )
                perc = 100
                if with_geom:
                    geom = geom_reference
            else:
                full = False
                full_total = False
                area = round(geom_intersection.area, 2)
                if with_geom:
                    geom = geom_intersection

            dict_observation["reference_features"][key_ref] = {
                "full": full,
                "area": area,
                "measure_type": feature_measure_type,
                feature_measure_type: round(intersection_measure, 2),
                "percentage": perc,
                "de9im": geometry.relate(geom_reference),
            }
            if version_date is not None:
                dict_observation["reference_features"][key_ref][VERSION_DATE] = (
                    version_date.strftime(DATE_FORMAT)
                )
            if with_geom:
                dict_observation["reference_features"][key_ref]["geometry"] = (
                    to_geojson(geom)
                )

        dict_observation["full"] = full_total
        if last_version_date is not None:
            dict_observation[LAST_VERSION_DATE] = last_version_date.strftime(
                DATE_FORMAT
            )
        if intersected:
            intersected_union = safe_unary_union(intersected)
        else:
            intersected_union = GeometryCollection()
        if measure_type == "area":
            geom_od = buffer_pos(
                buffer_neg(
                    safe_difference(geometry, intersected_union),
                    self.correction_distance,
                    mitre_limit=self.mitre_limit,
                ),
                self.correction_distance,
                mitre_limit=self.mitre_limit,
            )
        else:
            geom_od = safe_difference(geometry, intersected_union)
        if geom_od is not None:
            area_od = round(geom_od.area, 2)
            metric_od = round(_measure_value(geom_od, measure_type), 2)
            if area_od > 0 or metric_od > 0:
                dict_observation["reference_od"] = {"area": area_od}
                if measure_type != "area":
                    dict_observation["reference_od"][measure_type] = metric_od
                if with_geom:
                    dict_observation["reference_od"]["geometry"] = to_geojson(geom_od)
            if measure_type != "area":
                dict_observation["full"] = metric_od == 0 and bool(
                    dict_observation["reference_features"]
                )
        self.logger.feedback_debug(str(dict_observation))
        return dict_observation

    def get_difference_metrics_for_thematic_data(
        self,
        dict_processresults: Dict[InputId, Dict[float, Optional[ProcessResult]]] = None,
        thematic_data: AlignerFeatureCollection = None,
        diff_metric: DiffMetric = DiffMetric.SYMMETRICAL_AREA_CHANGE,
    ) -> Dict[InputId, Dict[float, float]]:
        """
        Calculates difference metrics for thematic elements across a series of distances.

        This method iterates through the processed results and compares the resulting
        geometries against their original thematic counterparts. It uses the
        spatial context of the reference data to calculate metrics like area
        change or symmetrical difference.

        Parameters
        ----------
        dict_processresults : Dict[InputId, Dict[float, Optional[ProcessResult]]], optional
            A nested dictionary where keys are thematic IDs and values are dictionaries
            mapping relevant distances to `ProcessResult`[] objects.
            Required if not provided via internal state.
        thematic_data : AlignerFeatureCollection, optional
            The collection containing the original thematic geometries.
            Defaults to `self.thematic_data`.
        diff_metric : DiffMetric, optional
            The metric used to quantify the geometric change.
            Defaults to `DiffMetric.SYMMETRICAL_AREA_CHANGE`.

        Returns
        -------
        Dict[InputId, Dict[float, float]]
            A nested dictionary where each thematic ID maps to a dictionary of
            distances and their corresponding calculated metric values.

        Raises
        ------
        ValueError
            If `dict_processresults` is not provided and cannot be resolved.

        Notes
        -----
        The calculation involves three primary geometric components:
        1. **Original Geometry**: The thematic feature before alignment.
        2. **Result Geometry**: The feature after alignment at a specific distance.
        3. **Reference Union**: The combined geometry of all reference data, used
           to contextualize the change.



        ```{mermaid}
        graph LR
            PR[Process Results] --> Calc[Metric Calculator]
            Orig[Original Geometries] --> Calc
            Ref[Reference Union] --> Calc
            Calc --> Output[Distance-Metric Map]
        ```

        Examples
        --------
        >>> metrics = aligner.get_difference_metrics_for_thematic_data(
        ...     dict_processresults=my_results,
        ...     diff_metric=DiffMetric.AREA_CHANGE
        ... )
        >>> # Get area change for feature 'A' at distance 0.5
        >>> print(metrics['A'][0.5])
        ```
        """
        if dict_processresults is None:
            raise ValueError("dict_processresults is required")
        if thematic_data is None:
            thematic_data = self.thematic_data
        diffs = {}
        for key, feat in thematic_data.features.items():
            original_geom = feat.geometry
            diffs[key] = get_geometry_difference_metrics_from_processresults(
                dict_processresult=dict_processresults[key],
                geom_thematic=original_geom,
                reference_union=self.reference_data.union,
                diff_metric=diff_metric,
            )

        return diffs

    def get_observation_properties(
        self,
        process_result,
    ):
        """
        Build observation properties for the current process result.

        Parameters
        ----------
        process_result : ProcessResult
            Processing output containing at least a `result` geometry.

        Returns
        -------
        dict
            Dictionary with derived observation flags, currently containing
            `FULL_ACTUAL_FIELD_NAME`.
        """
        return self.evaluator.get_observation_properties(
            aligner=self,
            process_result=process_result,
        )
