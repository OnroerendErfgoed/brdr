import json
from datetime import datetime
from typing import Dict, Optional, Any

import numpy as np

from brdr.aligner import Aligner
from brdr.be.grb.enums import GRBType
from brdr.be.grb.loader import GRBActualLoader
from brdr.be.grb.utils import get_affected_ids_by_grb_change
from brdr.configs import AlignerConfig
from brdr.constants import (
    LAST_VERSION_DATE,
    DATE_FORMAT,
    BASE_OBSERVATION_FIELD_NAME,
    DEFAULT_CRS,
    METADATA_FIELD_NAME,
    EVALUATION_FIELD_NAME,
)
from brdr.enums import FullReferenceStrategy, AlignerResultType, Evaluation
from brdr.loader import GeoJsonLoader
from brdr.logger import Logger


def update_featurecollection_to_actual_grb(
    featurecollection: Dict[str, Any],
    id_theme_fieldname: Optional[str] = None,
    base_metadata_field: str = METADATA_FIELD_NAME,
    grb_type: GRBType = GRBType.ADP,
    max_distance_for_actualisation: float = 2.0,
    max_predictions: int = -1,
    full_reference_strategy: FullReferenceStrategy = FullReferenceStrategy.NO_FULL_REFERENCE,
    multi_to_best_prediction: bool = True,
    feedback: Any = None,
    crs: str = DEFAULT_CRS,
    attributes: bool = True,
    max_workers: Optional[int] = None,
) -> Dict[str, Any]:
    """
    Update a thematic feature collection to the most recent version of the GRB.

    This function automates the alignment process by first identifying which
    features are 'affected' by changes in the GRB reference layer within a
    specific timeframe. It then calculates and evaluates new alignment
    predictions for those affected features.

    Parameters
    ----------
    featurecollection : Dict[str, Any]
        The thematic data as a GeoJSON-like dictionary.
    id_theme_fieldname : str, optional
        The property field name containing the unique ID for each feature.
    base_metadata_field : str, default OBSERVATION_FIELD_NAME
        Name of the attribute field that stores the existing alignment observation
        (JSON string). This is used to determine the last version date.
    grb_type : GRBType, default GRBType.ADP
        The type of GRB reference data to align against.
    max_distance_for_actualisation : float, default 2.0
        Maximum distance (in meters) to search for potential alignments.
        The function checks intervals of 0.1m up to this value.
    max_predictions : int, default -1
        Maximum number of alignment predictions to return. -1 returns all.
    full_reference_strategy : FullReferenceStrategy, default NO_FULL_REFERENCE
        Determines the prediction score when evaluating predictions, so predictions that are based on full-reference-geometries can be preferred.
    multi_to_best_prediction : bool, default True. Only useful when max_predictions=1
        If True, the prediction with the best prediction score is returned when multiple predictions are found
        If False, the original geometry is returned when multiple predictions are found
    feedback : Any, optional
        Feedback object (e.g., QgsFeedback) for progress reporting and logging (in QGIS).
    crs : str, default DEFAULT_CRS
        The Coordinate Reference System for processing.
    attributes : bool, default True
        Whether to include the original thematic attributes in the result.
    max_workers : int, optional
        The number of parallel threads to use for processing.

    Returns
    -------
    Dict[str, Any]
        A GeoJSON-like dictionary containing the evaluated predictions and
        updated alignment observations.

    Notes
    -----
    The function follows a specific lifecycle to optimize performance:

    1. **Initialization**: Loads thematic data and the actual GRB reference data.
    2. **Temporal Analysis**: Extracts the `last_version_date` from the features'
       observations to determine the relevant GRB change-window.
    3. **Spatial Filtering**: Uses `get_affected_by_grb_change` to isolate only
       those geometries where the underlying GRB has actually changed.
    4. **Alignment**: Executes the `Aligner.evaluate` logic only on 'affected'
       features.



    Raises
    ------
    ValueError
        If the CRS is unsupported or if the input feature collection is malformed.
    """
    last_version_date = None
    logger = Logger(feedback)

    # Initiate a Aligner to reference thematic features to the actual borders
    aligner_config = AlignerConfig()
    aligner_config.max_workers = max_workers
    aligner = Aligner(crs=crs, feedback=feedback, config=aligner_config)
    aligner.load_thematic_data(
        GeoJsonLoader(_input=featurecollection, id_property=id_theme_fieldname)
    )
    aligner.load_reference_data(
        GRBActualLoader(grb_type=grb_type, partition=1000, aligner=aligner)
    )
    rd_step = 10
    relevant_distances = [
        round(k, 1)
        for k in np.arange(
            0, max_distance_for_actualisation * 100 + rd_step, rd_step, dtype=int
        )
        / 100
    ]

    for id_theme, feature in aligner.thematic_data.features.items():
        try:
            if not base_metadata_field is None:
                base_observation_string = feature.properties[base_metadata_field]
                base_observation = json.loads(base_observation_string)

                logger.feedback_debug("observation: " + str(base_observation))
                try:
                    logger.feedback_debug(str(feature.properties))
                    if (
                        LAST_VERSION_DATE in base_observation
                        and base_observation[LAST_VERSION_DATE] is not None
                        and base_observation[LAST_VERSION_DATE] != ""
                    ):
                        str_lvd = base_observation[LAST_VERSION_DATE]
                        lvd = datetime.strptime(str_lvd, DATE_FORMAT).date()
                        if last_version_date is None or lvd < last_version_date:
                            last_version_date = lvd
                except Exception:
                    logger.feedback_debug(
                        f"Problem with {LAST_VERSION_DATE}. No brdr_observation (- json-attribute-field) loaded for id {id_theme}"
                    )
            else:
                logger.feedback_debug(
                    f"No brdr_observation (- json-attribute-field) loaded for id {str(id_theme)}"
                )
                last_version_date = None
        except:
            logger.feedback_debug(
                f"No brdr_observation (- json-attribute-field) loaded for id {str(id_theme)}"
            )
            last_version_date = None

    # als lastversiondate nog altijd 'now' is dan is er eigenlijk geen versiedate aanwezig in de data, en dan zetten we alle features op affected
    if last_version_date is not None:
        datetime_start = last_version_date
        datetime_end = datetime.now().date()
        thematic_geometries = {
            key: feat.geometry for key, feat in aligner.thematic_data.features.items()
        }

        affected = get_affected_ids_by_grb_change(
            thematic_geometries=thematic_geometries,
            grb_type=grb_type,
            date_start=datetime_start,
            date_end=datetime_end,
            one_by_one=False,
            geometry_thematic_union=aligner.thematic_data.union,
            border_distance=max_distance_for_actualisation,
            crs=aligner.crs,
        )
        logger.feedback_info(
            "Number of possible affected thematic geometries during timespan: "
            + str(len(affected))
        )
        if len(affected) == 0:
            logger.feedback_info(
                "No change detected in referencelayer during timespan. Script is finished"
            )
    else:
        affected = list(aligner.thematic_data.features.keys())

    aligner_result = aligner.process(
        thematic_ids=affected, relevant_distances=[max_distance_for_actualisation]
    )
    process_results = aligner_result.get_results(aligner=aligner)
    affected_and_changeable = []
    affected_and_not_changed = []
    for k, v in process_results.items():
        if not v[max_distance_for_actualisation]["result_diff"].is_empty:
            affected_and_changeable.append(k)
        else:
            affected_and_not_changed.append(k)

    # EXECUTE evaluation
    aligner_result = aligner.evaluate(
        relevant_distances=relevant_distances,
        thematic_ids=affected_and_changeable,
        metadata_field=BASE_OBSERVATION_FIELD_NAME,
        full_reference_strategy=full_reference_strategy,
        max_predictions=max_predictions,
        multi_to_best_prediction=multi_to_best_prediction,
    )
    for k, v in aligner_result.results.items():
        if k in affected_and_not_changed:
            for dist in v.keys():
                v[dist]["properties"][EVALUATION_FIELD_NAME] = Evaluation.NO_CHANGE

    return aligner_result.get_results_as_geojson(
        aligner=aligner,
        result_type=AlignerResultType.EVALUATED_PREDICTIONS,
        add_metadata=True,
        add_original_attributes=attributes,
    )
