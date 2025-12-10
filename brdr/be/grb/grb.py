import json
from datetime import datetime

import numpy as np
from shapely.geometry import shape, Polygon

from brdr.aligner import Aligner
from brdr.be.grb.enums import GRBType
from brdr.be.grb.loader import GRBActualLoader
from brdr.be.grb.utils import get_affected_by_grb_change
from brdr.constants import (
    LAST_VERSION_DATE,
    DATE_FORMAT,
    FORMULA_FIELD_NAME,
    BASE_FORMULA_FIELD_NAME,
)
from brdr.enums import FullReferenceStrategy, AlignerResultType
from brdr.loader import DictLoader
from brdr.logger import Logger


def update_to_actual_grb(
    featurecollection,
    id_theme_fieldname,
    base_formula_field=FORMULA_FIELD_NAME,
    grb_type=GRBType.ADP,
    max_distance_for_actualisation=2,
    max_predictions=-1,
    full_strategy=FullReferenceStrategy.NO_FULL_REFERENCE,
    multi_to_best_prediction=True,
    feedback=None,
    attributes=True,
):
    """
    Function to update a thematic featurecollection to the most actual version of GRB.
    Important to notice that the featurecollection needs a 'formula' for the base-alignment.

    :param featurecollection: Thematic featurecollection
    :param id_theme_fieldname: property-fieldname that states which property has to be used as unique ID
    :param base_formula_field: Name of the property-field that holds the original/base formula of the geometry, that has to be compared with the actual formula.
    :param max_distance_for_actualisation: Maximum relevant distance that is used to search and evaluate resulting geometries. All relevant distance between 0 and this max_distance are used to search, with a interval of 0.1m.
    :param feedback:  (default None): a QGIS feedback can be added to push all the logging to QGIS
    :param attributes: (boolean, default=True): States of all original attributes has to be added to the result
    :return: featurecollection
    """
    logger = Logger(feedback)

    # Load featurecollection into a shapely_dict
    dict_thematic = {}
    dict_thematic_props = {}

    last_version_date = None
    for feature in featurecollection["features"]:
        id_theme = feature["properties"][id_theme_fieldname]
        try:
            geom = shape(feature["geometry"])
        except:
            geom = Polygon()
        dict_thematic[id_theme] = geom
        dict_thematic_props[id_theme] = feature["properties"]
        try:
            if not base_formula_field is None:
                base_formula_string = feature["properties"][base_formula_field]
                dict_thematic_props[id_theme][
                    BASE_FORMULA_FIELD_NAME
                ] = base_formula_string
                base_formula = json.loads(base_formula_string)

                logger.feedback_debug("formula: " + str(base_formula))
                try:
                    logger.feedback_debug(str(dict_thematic_props[id_theme]))
                    if (
                        LAST_VERSION_DATE in base_formula
                        and base_formula[LAST_VERSION_DATE] is not None
                        and base_formula[LAST_VERSION_DATE] != ""
                    ):
                        str_lvd = base_formula[LAST_VERSION_DATE]
                        lvd = datetime.strptime(str_lvd, DATE_FORMAT).date()
                        if last_version_date is None or lvd < last_version_date:
                            last_version_date = lvd
                except Exception:
                    logger.feedback_info(
                        f"Problem with {LAST_VERSION_DATE}. No brdr_formula (- json-attribute-field) loaded for id {str(id_theme)}"
                    )
            else:
                logger.feedback_info(
                    f"No brdr_formula (- json-attribute-field) loaded for id {str(id_theme)}"
                )
                last_version_date = None
        except:
            logger.feedback_info(
                f"No brdr_formula (- json-attribute-field) loaded for id {str(id_theme)}"
            )
            last_version_date = None

    # als lastversiondate nog altijd 'now' is dan is er eigenlijk geen versiedate aanwezig in de data, en dan zetten we alle features op affected
    if last_version_date is not None:
        datetime_start = last_version_date
        datetime_end = datetime.now().date()
        base_aligner_result = Aligner(feedback=feedback)
        base_aligner_result.load_thematic_data(DictLoader(dict_thematic))
        base_aligner_result.name_thematic_id = id_theme_fieldname

        affected, unaffected = get_affected_by_grb_change(
            dict_thematic=base_aligner_result.dict_thematic,
            grb_type=grb_type,
            date_start=datetime_start,
            date_end=datetime_end,
            one_by_one=False,
            geometry_thematic_union=base_aligner_result.get_thematic_union(),
            border_distance=max_distance_for_actualisation,
            crs=base_aligner_result.CRS,
        )
        logger.feedback_info(
            "Number of possible affected OE-thematic during timespan: "
            + str(len(affected))
        )
        if len(affected) == 0:
            logger.feedback_info(
                "No change detected in referencelayer during timespan. Script is finished"
            )
    else:
        unaffected = []
        affected = list(dict_thematic.keys())

    # Initiate a Aligner to reference thematic features to the actual borders
    actual_aligner = Aligner(feedback=feedback, max_workers=None)
    actual_aligner.load_thematic_data(
        DictLoader(data_dict=dict_thematic, data_dict_properties=dict_thematic_props)
    )
    actual_aligner.load_reference_data(
        GRBActualLoader(grb_type=grb_type, partition=1000, aligner=actual_aligner)
    )
    rd_step = 10
    relevant_distances = [
        round(k, 1)
        for k in np.arange(
            0, max_distance_for_actualisation * 100 + rd_step, rd_step, dtype=int
        )
        / 100
    ]
    # EXECUTE evaluation
    aligner_result = actual_aligner.evaluate(
        ids_to_evaluate=affected,
        base_formula_field=BASE_FORMULA_FIELD_NAME,
        max_predictions=max_predictions,
        relevant_distances=relevant_distances,
        full_strategy=full_strategy,
        multi_to_best_prediction=multi_to_best_prediction,
    )

    return aligner_result.get_results_as_geojson(aligner=actual_aligner,result_type=AlignerResultType.EVALUATED_PREDICTIONS,
        formula=True,
        attributes=attributes,
    )
