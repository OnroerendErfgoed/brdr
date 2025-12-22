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
from typing import Any
from typing import Dict
from typing import Iterable
from typing import List
from typing import Optional
from typing import TYPE_CHECKING
from typing import Union

import numpy as np
from shapely import make_valid
from shapely import to_geojson
from shapely.geometry.base import BaseGeometry

from brdr import __version__
from brdr.configs import AlignerConfig
from brdr.configs import ProcessorConfig
from brdr.constants import AREA_CHANGE
from brdr.constants import AREA_PERCENTAGE_CHANGE
from brdr.constants import DATE_FORMAT
from brdr.constants import DEFAULT_CRS
from brdr.constants import DIFF_AREA_FIELD_NAME
from brdr.constants import DIFF_METRIC
from brdr.constants import DIFF_PERCENTAGE_FIELD_NAME
from brdr.constants import EQUAL_REFERENCE_FEATURES_FIELD_NAME
from brdr.constants import EVALUATION_FIELD_NAME
from brdr.constants import FORMULA_FIELD_NAME
from brdr.constants import FULL_ACTUAL_FIELD_NAME
from brdr.constants import FULL_BASE_FIELD_NAME
from brdr.constants import ID_THEME_FIELD_NAME
from brdr.constants import LAST_VERSION_DATE
from brdr.constants import LENGTH_CHANGE
from brdr.constants import LENGTH_PERCENTAGE_CHANGE
from brdr.constants import NR_CALCULATION_FIELD_NAME
from brdr.constants import OD_ALIKE_FIELD_NAME
from brdr.constants import PREDICTION_COUNT
from brdr.constants import PREDICTION_SCORE
from brdr.constants import RELEVANT_DISTANCE_DECIMALS
from brdr.constants import RELEVANT_DISTANCE_FIELD_NAME
from brdr.constants import REMARK_FIELD_NAME
from brdr.constants import STABILITY
from brdr.constants import SYMMETRICAL_AREA_CHANGE
from brdr.constants import SYMMETRICAL_AREA_PERCENTAGE_CHANGE
from brdr.constants import VERSION_DATE
from brdr.constants import ZERO_STREAK
from brdr.enums import AlignerResultType
from brdr.enums import DiffMetric
from brdr.enums import Evaluation
from brdr.enums import FullReferenceStrategy
from brdr.enums import ProcessRemark
from brdr.feature_data import AlignerFeatureCollection
from brdr.geometry_utils import buffer_neg
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
from brdr.typings import InputId
from brdr.utils import (
    coverage_ratio,
)
from brdr.utils import determine_stability
from brdr.utils import get_geojsons_from_process_results
from brdr.utils import get_geometry_difference_metrics_from_processresult
from brdr.utils import get_geometry_difference_metrics_from_processresults
from brdr.utils import is_brdr_formula
from brdr.utils import union_process_result
from brdr.utils import write_geojson

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
    results : Dict[ThematicId, Dict[float, Optional[ProcessResult]]]
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
        process_results : Dict[ThematicId, Dict[float, Optional[ProcessResult]]]
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
        Dict[ThematicId, Dict[float, Optional[ProcessResult]]]
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
                    (SYMMETRICAL_AREA_CHANGE, DIFF_METRIC.SYMMETRICAL_AREA_CHANGE),
                    (
                        SYMMETRICAL_AREA_PERCENTAGE_CHANGE,
                        DIFF_METRIC.SYMMETRICAL_AREA_PERCENTAGE_CHANGE,
                    ),
                    (AREA_CHANGE, DIFF_METRIC.AREA_CHANGE),
                    (AREA_PERCENTAGE_CHANGE, DIFF_METRIC.AREA_PERCENTAGE_CHANGE),
                    (LENGTH_CHANGE, DIFF_METRIC.LENGTH_CHANGE),
                    (LENGTH_PERCENTAGE_CHANGE, DIFF_METRIC.LENGTH_PERCENTAGE_CHANGE),
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
                    if res and EVALUATION_FIELD_NAME in res["properties"]
                }
                for tid, rd_dict in self.results.items()
            }
        else:
            raise ValueError(f"Unknown result type: {result_type}")

    def get_results_as_geojson(
        self,
        aligner: "Aligner",
        result_type: AlignerResultType = AlignerResultType.PROCESSRESULTS,
        formula: bool = False,
        attributes: bool = False,
    ) -> Dict[str, Any]:
        """
        Converts the results into a GeoJSON FeatureCollection format.

        Parameters
        ----------
        aligner : Aligner
            The 'Aligner' object providing CRS and thematic metadata.
        result_type : AlignerResultType
            The type of results to export. Defaults to PROCESSRESULTS.
        formula : bool
            If True, includes the descriptive comparison formula in properties.
        attributes : bool
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

                if attributes and feature and feature.properties:
                    prop_dictionary[theme_id][relevant_distance].update(
                        feature.properties
                    )

                if formula:
                    result_geom = process_result["result"]
                    formula_result = process_result.get(
                        "formula"
                    ) or aligner.compare_to_reference(result_geom)
                    prop_dictionary[theme_id][relevant_distance][FORMULA_FIELD_NAME] = (
                        json.dumps(formula_result)
                    )

        return get_geojsons_from_process_results(
            results,
            crs=aligner.crs,
            id_field=aligner.thematic_data.id_fieldname,
            series_prop_dict=prop_dictionary,
        )

    def save_results(
        self,
        aligner: "Aligner",
        path: str,
        result_type: AlignerResultType = AlignerResultType.PROCESSRESULTS,
        formula: bool = False,
        attributes: bool = False,
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
        formula : bool
            Whether to include alignment formulas in the output.
        attributes : bool
            Whether to include original feature attributes.

        Notes
        -----
        The output files follow the naming convention: `{result_type}_{name}.geojson`.
        """
        fcs = self.get_results_as_geojson(
            aligner=aligner,
            result_type=result_type,
            formula=formula,
            attributes=attributes,
        )
        for name, fc in fcs.items():
            file_name = f"{result_type.value}_{name}.geojson"
            write_geojson(os.path.join(path, file_name), fc)


def _get_observations_from_formula(processResult: ProcessResult) -> List[Dict]:
    formula = processResult["formula"]
    actuation_metadata = processResult["metadata"]["actuation"]
    observation_time = datetime.now().strftime("%Y-%m-%dT%H:%M:%S%z")
    sensor_uuid = uuid.uuid4().hex
    result_id = actuation_metadata["result"]
    observation_metadata = {
        "type": "sosa:Observation",
        "has_feature_of_interest": result_id,
        "made_by_sensor": f"brdrid:sensors/{sensor_uuid}",
        "result_time": observation_time,
    }

    def get_reference_from_actuation(ref_id):
        for ref in actuation_metadata["reference_geometries"]:
            if ref["id"] == ref_id:
                return ref
        raise ValueError("reference geometry not found in actuation metadata")

    observations = []
    for ref_id, observations_dict in formula["reference_features"].items():
        reference = get_reference_from_actuation(ref_id)
        if area:= observations_dict.get("area"):
            observations.append(
                {
                    **observation_metadata,
                    "id": f"brdrid/observations/{uuid.uuid4().hex}",
                    "observed_property": "brdr:area_overlap",
                    "result": {"value": area, "type": "float"},
                    "used_procedure": "brdr:observation_procedure_area_overlap",
                    "used": reference,
                }
            )
        if percentage:= observations_dict.get("percentage"):
            observations.append(
                {
                    **observation_metadata,
                    "id": f"brdrid/observations/{uuid.uuid4().hex}",
                    "observed_property": "brdr:area_overlap_percentage",
                    "result": {"value": percentage, "type": "float"},
                    "used_procedure": "brdr:observation_procedure_area_overlap_percentage",
                    "used": reference,
                }
            )
        if full:= observations_dict.get("full"):
            observations.append(
                {
                    **observation_metadata,
                    "id": f"brdrid/observations/{uuid.uuid4().hex}",
                    "observed_property": "brdr:area_overlap_full",
                    "result": {"value": full, "type": "boolean"},
                    "used_procedure": "brdr:observation_procedure_area_overlap_full",
                    "used": reference,
                }
            )

    if area:= formula.get("area"):
        observations.append(
            {
                **observation_metadata,
                "id": f"brdrid/observations/{uuid.uuid4().hex}",
                "observed_property": "brdr:area",
                "result": {"value": area, "type": "float"},
                "used_procedure": "brdr:observation_procedure_area",
                "used": {
                    "id": result_id,
                    "type": f"geo:Geometry",
                    "version_date": observation_time,
                },
            }
        )
    if area_od:= formula.get("area_od", {}).get("area"):
        observations.append(
            {
                **observation_metadata,
                "id": f"brdrid/observations/{uuid.uuid4().hex}",
                "observed_property": "brdr:area_od",
                "result": {"value": area_od, "type": "float"},
                "used_procedure": "brdr:observation_procedure_area_od",
                "used": {
                    "id": result_id,
                    "type": f"geo:Geometry",
                    "version_date": observation_time,
                },
            }
        )
    if full:= formula.get("full"):
        observations.append(
            {
                **observation_metadata,
                "id": f"brdrid/observations/{uuid.uuid4().hex}",
                "observed_property": "brdr:area_overlap_full",
                "result": {"value": full, "type": "bool"},
                "used_procedure": "brdr:observation_procedure_area_overlap_full",
            }
        )
    assert len(observations) > 0


def aligner_metadata_decorator(f):
    def inner_func(aligner, *args, **kwargs):
        assert isinstance(aligner, Aligner)
        response: AlignerResult = f(aligner, *args, **kwargs)
        if aligner.add_formula:
            for thematic_id, rd_res in response.results.items():
                for rd, res in rd_res.items():
                    if res["result"]:
                        formula = aligner.compare_to_reference(res["result"])
                        res["formula"] = formula
        if aligner.log_metadata:
            # generate uuid for actuation
            actuation_id = "brdrid:actuations/" + uuid.uuid4().hex
            processor_id = aligner.processor.processor_id.value
            processor_name = type(aligner.processor).__name__
            reference_data = aligner.reference_data
            reference_features = reference_data.features.values()
            reference_geometries = [
                {
                    "id": feature.data_id,
                    "type": f"geo:{feature.geometry.geom_type}",
                    "version_date": reference_data.source.get(VERSION_DATE, ""),
                    "identifier": {
                        "id": feature.brdr_id,
                        "type": "adms:Identifier",
                    },
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
                result_hash = hashlib.sha1(result["result"].wkt.encode()).hexdigest()
                result["metadata"] = {}
                result["metadata"]["actuation"] = {
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
                if result["formula"]:
                    result["metadata"]["observations"] = _get_observations_from_formula(
                        result
                    )

        return response

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
    add_formula: bool
        If True, process result formulas will be computed by default.
    processor : BaseProcessor or AlignerGeometryProcessor
        The geometric processor used for alignment calculations.
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
    >>> aligner.load_thematic_data(my_features)
    >>> results = aligner.align(relevant_distances=[0.5, 1.0])
    """

    def __init__(
        self,
        *,
        processor: Optional[BaseProcessor] = None,
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
        self.add_formula = config.add_formula
        self.processor = (
            processor
            if processor
            else AlignerGeometryProcessor(ProcessorConfig(), feedback)
        )
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

    @aligner_metadata_decorator
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
        thematic_ids : List[ThematicId], optional
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

        def process_geom_for_rd(geometry, relevant_distance):
            return self.processor.process(
                correction_distance=self.correction_distance,
                reference_data=self.reference_data,
                mitre_limit=self.mitre_limit,
                input_geometry=geometry,
                relevant_distance=relevant_distance,
                thematic_data=self.thematic_data,
            )

        def run_process(executor: ThreadPoolExecutor = None):
            for thematic_id in thematic_ids:
                geom = self.thematic_data.features.get(thematic_id).geometry
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
        relevant_distances: Optional[Union[List[float], np.ndarray]] = None,
        *,
        thematic_ids: Optional[List[InputId]] = None,
        diff_metric: Optional[DIFF_METRIC] = None,
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
        thematic_ids : List[ThematicId], optional
            Specific thematic IDs to process. If None, all loaded thematic
            geometries are used.
        diff_metric : DIFF_METRIC, optional
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

        if thematic_ids is None:
            thematic_ids = self.thematic_data.features.keys()
        if any(
            id_to_predict not in self.thematic_data.features.keys()
            for id_to_predict in thematic_ids
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
        process_all_at_once = True
        if process_all_at_once:
            aligner_result = self.process(
                thematic_ids=thematic_ids,
                relevant_distances=rd_prediction,
            )
            process_results = aligner_result.results

        # else:
        #     process_results = {}
        #     futures={}
        #
        #     def _predict_geom(theme_id,geom,_rd_prediction):
        #
        #         def _check_interval_stability(
        #             res_start: ProcessResult, res_end: ProcessResult
        #         ) -> bool:
        #
        #             return geometric_equality(
        #                 res_start["result"],
        #                 res_end["result"],
        #                 correction_distance=self.correction_distance,
        #                 mitre_limit=self.mitre_limit,
        #             )
        #
        #         def _process_result(_theme_id, _relevant_distances):
        #             aligner_result = self.process(
        #                 thematic_ids=[_theme_id],
        #                 relevant_distances=_relevant_distances,
        #                 max_workers=self.max_workers)
        #             return aligner_result.results[_theme_id][_relevant_distances[0]]
        #
        #         def _process_wrapper(x: float):
        #             return _process_result(_theme_id=theme_id, _relevant_distances=[x])
        #
        #         non_stable_points, cache = recursive_stepwise_interval_check(
        #             f_value=_process_wrapper,
        #             f_condition=_check_interval_stability,
        #             discrete_list=_rd_prediction,
        #             initial_sample_size=3,
        #         )
        #         interpolated_cache = create_full_interpolated_dataset(_rd_prediction, cache)
        #         return interpolated_cache
        #
        #     def run_prediction(executor: ThreadPoolExecutor = None):
        #
        #         for key in thematic_ids:
        #             geom=self.thematic_data.features.get(key).geometry
        #             self.logger.feedback_info(
        #                 f"thematic id {str(key)} predictions calculated with "
        #                 f"relevant distances (m) [{str(relevant_distances)}]"
        #             )
        #             process_results[key] = {}
        #             try:
        #                 fn = _predict_geom
        #                 if executor:
        #                     futures[key] = executor.submit(fn, key,geom,rd_prediction)
        #                 else:
        #                     process_results[key] = fn(key,geom,rd_prediction)
        #             except ValueError as e:
        #                 self.logger.feedback_warning(
        #                     f"error for thematic id {str(key)} processed with "
        #                     f"relevant distances (m) [{str(relevant_distances)}]"
        #                 )
        #                 process_results[key] = None
        #                 self.logger.feedback_warning(str(e))
        #
        #     if self.max_workers == -1:
        #         run_prediction()
        #     else:
        #         with ThreadPoolExecutor(self.max_workers) as prediction_executor:
        #             run_prediction(prediction_executor)
        #             self.logger.feedback_debug("waiting all started thematic calculations")
        #             wait(list(futures.values()))
        #             for thematic_id, future in futures.items():
        #                 process_results[thematic_id] = future.result()

        if diff_metric is None:
            diff_metric = self.diff_metric
        diffs_dict = {}
        for theme_id, process_result in process_results.items():
            diffs = get_geometry_difference_metrics_from_processresults(
                process_result,
                self.thematic_data.features.get(theme_id).geometry,
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
        relevant_distances: Optional[Iterable[float]] = None,
        *,
        thematic_ids: Optional[List[InputId]] = None,
        base_formula_field: str = FORMULA_FIELD_NAME,
        full_reference_strategy: FullReferenceStrategy = FullReferenceStrategy.NO_FULL_REFERENCE,
        max_predictions: int = -1,
        multi_to_best_prediction: bool = True,
    ) -> AlignerResult:
        """
        Evaluates input geometries against predictions using formula matching
        and selection strategies.

        This method identifies the best candidate geometries for alignment by
        comparing predicted geometries against the original attributes and
        reference data. It tags results with evaluation labels (e.g.,
        PREDICTION_UNIQUE) to facilitate decision-making.

        Parameters
        ----------
        relevant_distances : Iterable[float], optional
            Distances to evaluate. Defaults to 0.0m to 3.0m (step 0.1m).
        thematic_ids : List[ThematicId], optional
            List of IDs to evaluate. If None, all loaded thematic features
            are processed. Features not in this list are marked as NOT_EVALUATED.
        base_formula_field : str, optional
            The field name containing the base comparison formula.
            Defaults to FORMULA_FIELD_NAME.
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
        1. **Formula Match**: If a prediction matches the base formula exactly,
           it is prioritized (score 100).
        2. **Fullness**: If a strategy is set, 'full' predictions get a score boost.
        3. **Scoring**: Candidates are ranked by their prediction score.
        4. **Selection**: Based on `max_predictions`, the best candidates are selected.



        ```{mermaid}
        graph TD
            Start[Start Evaluation] --> Pred[Generate Predictions]
            Pred --> Match{Formula Match?}
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
        >>> results = aligner.evaluate(
        ...     thematic_ids=["id_01"],
        ...     max_predictions=1,
        ...     full_reference_strategy=FullReferenceStrategy.ONLY_FULL_REFERENCE
        ... )
        """

        if thematic_ids is None:
            thematic_ids = self.thematic_data.features.keys()
        if any(
            id_to_evaluate not in self.thematic_data.features.keys()
            for id_to_evaluate in thematic_ids
        ):
            raise ValueError("not all ids are found in the thematic data")
        if relevant_distances is None:
            relevant_distances = [
                round(k, RELEVANT_DISTANCE_DECIMALS)
                for k in np.arange(0, 310, 10, dtype=int) / 100
            ]

        aligner_result = self.predict(
            thematic_ids=thematic_ids,
            relevant_distances=relevant_distances,
            diff_metric=self.diff_metric,
        )
        process_results_evaluated = aligner_result.get_results(aligner=self)
        process_results_evaluated_predictions = aligner_result.get_results(
            aligner=self, result_type=AlignerResultType.PREDICTIONS
        )
        process_results_evaluated = deepcopy(process_results_evaluated)

        for theme_id, feat in self.thematic_data.features.items():
            original_geometry = feat.geometry
            # Features are split up in 2 groups: TO_EVALUATE and NOT_TO_EVALUATE (original returned)
            # The evaluated features will be split up:
            #   *No prediction available
            #   *Predictions available

            # PART 1: TO_EVALUATE
            if theme_id in thematic_ids:
                if theme_id not in process_results_evaluated_predictions.keys():
                    # No predictions available
                    relevant_distance = round(0, RELEVANT_DISTANCE_DECIMALS)
                    props = self._evaluate(
                        id_theme=theme_id,
                        geom_predicted=original_geometry,
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
                    process_results_evaluated[theme_id][relevant_distance] = (
                        union_process_result(
                            {
                                "result": original_geometry,
                                "properties": props,
                            }
                        )
                    )
                    continue

                # When there are predictions available
                dict_predictions_results = process_results_evaluated_predictions[
                    theme_id
                ]
                scores = []
                distances = []
                predictions = []
                formula_match = False
                for dist in sorted(dict_predictions_results.keys()):
                    props = self._evaluate(
                        id_theme=theme_id,
                        geom_predicted=dict_predictions_results[dist]["result"],
                        base_formula_field=base_formula_field,
                        formula=dict_predictions_results[dist].get("formula"),
                    )
                    props.update(
                        process_results_evaluated_predictions[theme_id][dist][
                            "properties"
                        ]
                    )

                    full = props[FULL_ACTUAL_FIELD_NAME]
                    if (
                        full_reference_strategy
                        == FullReferenceStrategy.ONLY_FULL_REFERENCE
                        and not full
                    ):
                        # this prediction is ignored
                        continue
                    if (
                        props[EVALUATION_FIELD_NAME]
                        == Evaluation.TO_CHECK_NO_PREDICTION
                        and props[PREDICTION_COUNT] == 1
                    ):
                        props[EVALUATION_FIELD_NAME] = Evaluation.PREDICTION_UNIQUE
                    elif (
                        props[EVALUATION_FIELD_NAME]
                        == Evaluation.TO_CHECK_NO_PREDICTION
                        and props[PREDICTION_COUNT] > 1
                    ):
                        props[EVALUATION_FIELD_NAME] = (
                            Evaluation.TO_CHECK_PREDICTION_MULTI
                        )
                    elif (
                        props[EVALUATION_FIELD_NAME]
                        != Evaluation.TO_CHECK_NO_PREDICTION
                    ):
                        # this prediction has a equality based on formula so the rest is not checked anymore
                        formula_match = True
                        props[PREDICTION_SCORE] = 100
                        scores = []
                        distances = []
                        predictions = []
                        scores.append(props[PREDICTION_SCORE])
                        distances.append(dist)
                        process_results_evaluated_predictions[theme_id][dist][
                            "properties"
                        ] = props
                        predictions.append(
                            process_results_evaluated_predictions[theme_id][dist]
                        )
                        continue
                    if full:
                        if (
                            full_reference_strategy
                            != FullReferenceStrategy.NO_FULL_REFERENCE
                        ):
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
                    process_results_evaluated_predictions[theme_id][dist][
                        "properties"
                    ] = props
                    predictions.append(
                        process_results_evaluated_predictions[theme_id][dist]
                    )

                # get max amount of best-scoring predictions
                best_ix = sorted(
                    range(len(scores)), reverse=True, key=lambda i: scores[i]
                )
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
                        remarks.append(
                            ProcessRemark.MULTIPLE_PREDICTIONS_ORIGINAL_RETURNED
                        )
                        props[REMARK_FIELD_NAME] = remarks
                        process_results_evaluated[theme_id][relevant_distance] = (
                            union_process_result(
                                {
                                    "result": original_geometry,
                                    "properties": props,
                                }
                            )
                        )
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
                        geom_predicted=original_geometry,
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
                    process_results_evaluated[theme_id][relevant_distance] = (
                        union_process_result(
                            {
                                "result": original_geometry,
                                "properties": props,
                            }
                        )
                    )
            else:
                # PART 2: NOT EVALUATED
                relevant_distance = round(0, RELEVANT_DISTANCE_DECIMALS)
                process_results_evaluated[theme_id] = {}
                props = self._evaluate(
                    id_theme=theme_id,
                    geom_predicted=original_geometry,
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
                process_results_evaluated[theme_id][relevant_distance] = (
                    union_process_result(
                        {
                            "result": original_geometry,
                            "properties": props,
                        }
                    )
                )
        return AlignerResult(process_results_evaluated)

    def compare_to_reference(
        self, geometry: BaseGeometry, with_geom: bool = False
    ) -> Dict[str, Any]:
        """
        Calculates formula-related information based on the input geometry.

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
            Stats --> Dict[Final Formula Dictionary]
            OD --> Dict
        ```

        Examples
        --------
        >>> info = aligner.compare_to_reference(my_geom, with_geom=True)
        >>> print(info["full"])
        >>> print(info["reference_features"].keys())  # IDs of intersected refs
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

    def get_difference_metrics_for_thematic_data(
        self,
        dict_processresults: Dict[
            InputId, Dict[float, Optional[ProcessResult]]
        ] = None,
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
        dict_processresults : Dict[ThematicId, Dict[float, Optional[ProcessResult]]], optional
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
        Dict[ThematicId, Dict[float, float]]
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

    def _evaluate(
        self,
        id_theme,
        geom_predicted,
        base_formula_field=FORMULA_FIELD_NAME,
        formula=None,
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
        actual_formula = formula or self.compare_to_reference(geom_predicted)
        if actual_formula is None or geom_predicted is None or geom_predicted.is_empty:
            properties[EVALUATION_FIELD_NAME] = Evaluation.TO_CHECK_NO_PREDICTION
            properties[FULL_ACTUAL_FIELD_NAME] = False
            return properties
        properties[FULL_ACTUAL_FIELD_NAME] = actual_formula["full"]
        properties[FORMULA_FIELD_NAME] = json.dumps(actual_formula)

        try:
            base_formula = json.loads(
                self.thematic_data.features.get(id_theme).properties[base_formula_field]
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
