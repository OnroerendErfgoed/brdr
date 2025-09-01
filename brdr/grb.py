import json
import logging
from copy import deepcopy
from datetime import date
from datetime import datetime

import numpy as np
from shapely import intersects
from shapely.geometry import shape, Polygon

from brdr.aligner import Aligner
from brdr.constants import (
    DEFAULT_CRS,
    LAST_VERSION_DATE,
    DATE_FORMAT,
    VERSION_DATE,
    FORMULA_FIELD_NAME,
    BASE_FORMULA_FIELD_NAME,
    GRB_GENERIC_ID,
)
from brdr.constants import DOWNLOAD_LIMIT
from brdr.constants import GRB_FEATURE_URL
from brdr.constants import GRB_FISCAL_PARCELS_URL
from brdr.constants import GRB_MAX_REFERENCE_BUFFER
from brdr.constants import GRB_PARCEL_ID
from brdr.constants import GRB_VERSION_DATE
from brdr.enums import GRBType, AlignerResultType, FullStrategy
from brdr.geometry_utils import buffer_pos, safe_intersection, safe_unary_union
from brdr.geometry_utils import create_donut
from brdr.geometry_utils import features_by_geometric_operation
from brdr.geometry_utils import get_bbox
from brdr.loader import GeoJsonLoader, DictLoader
from brdr.logger import Logger
from brdr.utils import geojson_to_dicts
from brdr.utils import get_collection
from brdr.utils import get_collection_by_partition

log = logging.getLogger(__name__)

datetime_format_TZ = "%Y-%m-%dT%H:%M:%SZ"


def is_grb_changed(
    geometry,
    grb_type=GRBType.ADP,
    date_start=date.today(),
    date_end=date.today(),
    border_distance=0,
    crs=DEFAULT_CRS,
):
    """
       checks if a geometry is possibly affected by changes in the reference layer
       during a specified timespan
    Args:
        geometry: Geometry to check on GRB-changes
        grb_type:  Type of GRB (parcels, buildings,artwork,...) to check
        date_start: Start of timespan to check if GRB-changes has occurred
        date_end: End of timespan to check if GRB-changes has occurred
        border_distance: Distance that can be used to only check the 'border' of
            the geometry, so 'big' geometries with internal parcel-updates are not
             affected
        (Default:0, indicating that the full geometry is checked fot GRB-changes)
        crs: Coordinate reference system to use

    Returns: Boolean, indicating if GRB is changed underneath the geometry

    """

    last_version_date = get_last_version_date(
        geometry, grb_type=grb_type, crs=crs, border_distance=border_distance
    )
    if last_version_date is None:
        return None
    if date_start <= last_version_date <= date_end:
        return True
    return False


def get_affected_by_grb_change(
    dict_thematic,
    grb_type=GRBType.ADP,
    date_start=date.today(),
    date_end=date.today(),
    one_by_one=False,
    border_distance=0,
    geometry_thematic_union=None,
    crs=DEFAULT_CRS,
):
    """
    Get a list of affected and unaffected IDs by GRB-changes in a
    specific timespan

    Args:
        dict_thematic: dictionary if thematicID & Geometry
        grb_type: Type of GRB: parcels, buildings,...
        date_start: start-date to check changes in GRB
        date_end: end-date to check changes in GRB
        one_by_one: parameter to choose the methodology to check changes:
            * True: Every thematic geometry is checked individually (loop)
            * False: All GRB-parcels intersecting the thematic dictionary is checked
                at once
        border_distance: Distance that can be used to only check the 'border' of  the
            geometry, so 'big' geometries with internal parcel-updates are not affected
            (Default:0, indicating that the full geometry is checked fot GRB-changes)
    Returns:
        dictionary of affected geometries

    """

    affected = []
    unaffected = []
    if one_by_one:
        for key in dict_thematic:
            geom = dict_thematic[key]
            if is_grb_changed(
                geom, grb_type, date_start, date_end, border_distance=border_distance
            ):
                affected.append(key)
            else:
                unaffected.append(key)
        return affected, unaffected
    else:
        # Temporal filter on VERDATUM
        if geometry_thematic_union is None:
            geometry_thematic_union = safe_unary_union(list(dict_thematic.values()))
        coll_changed_grb, name_reference_id = get_collection_grb_actual(
            buffer_pos(geometry_thematic_union, GRB_MAX_REFERENCE_BUFFER),
            grb_type=grb_type,
            partition=1000,
            date_start=date_start,
            date_end=date_end,
            crs=crs,
        )
        dict_changed_grb, dict_changed_grb_properties = geojson_to_dicts(
            coll_changed_grb, name_reference_id
        )
        # if border_distance>0:
        geom_to_check = buffer_pos(
            create_donut(geometry_thematic_union, border_distance),
            GRB_MAX_REFERENCE_BUFFER,
        )
        grb_intersections = features_by_geometric_operation(
            list(dict_changed_grb.values()),
            list(dict_changed_grb.keys()),
            [geom_to_check],
            predicate="intersects",
        )
        dict_changed_grb = {key: dict_changed_grb[key] for key in grb_intersections}

        if len(dict_changed_grb) == 0:
            logging.info(
                f"No detected changes for thematic geometry in timespan (border distance: {str(border_distance)})"
            )
            return affected, list(dict_thematic.keys())  # empty affected dict
        logging.info(
            f"Changed parcels in timespan with border_distance {str(border_distance)}: {str(len(dict_changed_grb))}"
        )
        thematic_intersections = features_by_geometric_operation(
            list(dict_thematic.values()),
            list(dict_thematic.keys()),
            list(dict_changed_grb.values()),
            predicate="intersects",
        )
        logging.info("Number of filtered features: " + str(len(thematic_intersections)))
        for key, geom in dict_thematic.items():
            (
                affected.append(key)
                if key in thematic_intersections
                else unaffected.append(key)
            )
    return affected, unaffected


def get_last_version_date(
    geometry,
    grb_type=GRBType.ADP,
    crs=DEFAULT_CRS,
    limit=DOWNLOAD_LIMIT,
    border_distance=0,
):
    """
    Retrieves the date of the last version for a specific geographic area within
    GRB (parcels, buildings,...).

    This function queries the GRB-API to find the most recent version-date (=last
    update  of object) for reference data of the specified `grb_type` (e.g., ADP,
    GBG, KNW) within the boundary of the provided `geometry`.

    Args:
    geometry (BaseGeometry): A Shapely geometry representing the area of interest.
    grb_type (GRBType, optional): The type of GRB data to consider.
        Defaults to GRBType.ADP (administrative plots).

    Returns:
    str: The date of the last version for the specified GRB data type within the area,
         formatted as a string according to the GRB API response (usually YYYY-MM-DD).

    None: If no data is found for the given geometry and GRB type combination.
    """
    if border_distance > 0:
        geometry = create_donut(geometry, border_distance)
    bbox = get_bbox(geometry)
    actual_url = (
        GRB_FEATURE_URL + "/" + grb_type.name + "/items?"
        "limit=" + str(limit) + "&crs=" + crs + "&bbox-crs=" + crs + "&bbox=" + bbox
    )
    update_dates = []
    collection = get_collection(actual_url, limit)
    if "features" not in collection:
        return None
    for c in collection["features"]:
        if intersects(geometry, shape(c["geometry"])):
            versiondate = datetime.strptime(
                c["properties"][GRB_VERSION_DATE], DATE_FORMAT
            ).date()
            update_dates.append(versiondate)

    update_dates = sorted(update_dates, reverse=True)
    if len(update_dates) > 0:
        return update_dates[0]
    return None


def get_collection_grb_actual(
    geometry,
    grb_type=GRBType.ADP,
    partition=1000,
    limit=DOWNLOAD_LIMIT,
    crs=DEFAULT_CRS,
    date_start=None,
    date_end=None,
):
    url = (
        GRB_FEATURE_URL
        + "/"
        + grb_type.name
        + "/items?limit="
        + str(limit)
        + "&crs="
        + crs
    )
    if grb_type == GRBType.ADP:
        name_reference_id = GRB_PARCEL_ID
    else:
        name_reference_id = GRB_GENERIC_ID

    versiondate_filter = ""
    versiondate_filter_start = ""
    versiondate_filter_end = ""
    if date_start is not None:
        versiondate_filter_start = (
            GRB_VERSION_DATE + ">" + date_start.strftime(DATE_FORMAT)
        )
        versiondate_filter = versiondate_filter_start
    if date_end is not None:
        versiondate_filter_end = GRB_VERSION_DATE + "<" + date_end.strftime(DATE_FORMAT)
        versiondate_filter = versiondate_filter_end
    if not (date_start is None and date_end is None):
        versiondate_filter = versiondate_filter_start + " AND " + versiondate_filter_end
    if versiondate_filter != "":
        url = url + "&filter=" + versiondate_filter + "&filter-lang=cql-text"

    collection = get_collection_by_partition(
        url, geometry=geometry, partition=partition, limit=limit, crs=crs
    )
    return collection, name_reference_id


def get_collection_grb_fiscal_parcels(
    geometry,
    year=str(datetime.now().year),
    partition=1000,
    limit=DOWNLOAD_LIMIT,
    crs=DEFAULT_CRS,
):
    url = (
        GRB_FISCAL_PARCELS_URL + "/Adpf" + year + "/items?"
        "limit=" + str(limit) + "&crs=" + crs
    )
    return get_collection_by_partition(
        url, geometry=geometry, partition=partition, limit=limit, crs=crs
    )


def get_collection_grb_parcels_by_date(
    geometry,
    date,
    partition=1000,
    limit=DOWNLOAD_LIMIT,
    crs=DEFAULT_CRS,
):
    collection_year_after = get_collection_grb_fiscal_parcels(
        year=str(date.year),
        geometry=geometry,
        partition=partition,
        limit=limit,
        crs=crs,
    )
    # Filter on specific date: delete all features > specific_date
    # TODO: experimental loader; unclear if we have to use "year-1 & year" OR if we have to use "year & year + 1"
    collection_year_after_filtered = deepcopy(collection_year_after)
    logging.debug(len(collection_year_after_filtered["features"]))
    if (
        "features" in collection_year_after_filtered
        and len(collection_year_after_filtered["features"]) > 0
    ):
        removed_features = []
        for feature in collection_year_after_filtered["features"]:
            versiondate = datetime.strptime(
                feature["properties"][GRB_VERSION_DATE][:10], DATE_FORMAT
            ).date()
            if versiondate > date:
                removed_features.append(feature)
                collection_year_after_filtered["features"].remove(feature)
    logging.debug(len(collection_year_after_filtered["features"]))
    # if no features are removed, return the full collection of year_after
    if len(removed_features) == 0:
        return collection_year_after
    # if  features are removed, search for the features in year before
    collection_year_before = get_collection_grb_fiscal_parcels(
        year=str(date.year - 1),
        geometry=geometry,
        partition=partition,
        limit=limit,
        crs=crs,
    )
    kept_features = []
    if "features" in collection_year_before and len(collection_year_before) > 0:
        for feature in collection_year_before["features"]:
            for rf in removed_features:
                geom_feature = shape(feature["geometry"])
                geom_removed_feature = shape(rf["geometry"])
                if intersects(geom_feature, geom_removed_feature):
                    intersection = safe_intersection(geom_feature, geom_removed_feature)
                    if intersection.area > 1:
                        if feature not in kept_features:
                            kept_features.append(feature)

        # search for intersection and check if it more than x%
        # keep these features

        # add them to

    collection_specific_date = deepcopy(collection_year_after_filtered)
    filtered_features = collection_year_after_filtered["features"]
    specific_date_features = filtered_features + kept_features
    logging.debug(len(specific_date_features))
    collection_specific_date["features"] = specific_date_features

    return collection_specific_date


def update_to_actual_grb(
    featurecollection,
    id_theme_fieldname,
    base_formula_field=FORMULA_FIELD_NAME,
    grb_type=GRBType.ADP,
    max_distance_for_actualisation=2,
    max_predictions=-1,
    full_strategy=FullStrategy.NO_FULL,
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
    actual_aligner.evaluate(
        ids_to_evaluate=affected,
        base_formula_field=BASE_FORMULA_FIELD_NAME,
        max_predictions=max_predictions,
        relevant_distances=relevant_distances,
        full_strategy=full_strategy,
        multi_to_best_prediction=multi_to_best_prediction,
    )

    return actual_aligner.get_results_as_geojson(
        resulttype=AlignerResultType.EVALUATED_PREDICTIONS,
        formula=True,
        attributes=attributes,
    )


class GRBActualLoader(GeoJsonLoader):
    def __init__(self, grb_type: GRBType, aligner, partition: int = 1000):
        super().__init__()
        self.aligner = aligner
        self.grb_type = grb_type
        self.part = partition
        self.data_dict_source["source"] = grb_type.value
        self.versiondate_info = {"name": GRB_VERSION_DATE, "format": DATE_FORMAT}

    def load_data(self):
        if not self.aligner.dict_thematic:
            raise ValueError("Thematic data not loaded")
        geom_union = buffer_pos(
            self.aligner.get_thematic_union(), GRB_MAX_REFERENCE_BUFFER
        )
        collection, id_property = get_collection_grb_actual(
            grb_type=self.grb_type,
            geometry=geom_union,
            partition=self.part,
            crs=self.aligner.CRS,
        )
        self.id_property = id_property
        self.input = dict(collection)
        self.data_dict_source[VERSION_DATE] = datetime.now().strftime(DATE_FORMAT)
        self.aligner.logger.feedback_info(f"GRB downloaded: {self.grb_type}")
        return super().load_data()


class GRBFiscalParcelLoader(GeoJsonLoader):
    def __init__(self, year: str, aligner, partition=1000):
        super().__init__(_input=None, id_property=GRB_PARCEL_ID)
        self.aligner = aligner
        self.year = year
        self.part = partition
        self.data_dict_source["source"] = "Adpf"
        self.data_dict_source[VERSION_DATE] = datetime(int(year), 1, 1).strftime(
            DATE_FORMAT
        )
        self.versiondate_info = {"name": GRB_VERSION_DATE, "format": datetime_format_TZ}

    def load_data(self):
        if not self.aligner.dict_thematic:
            raise ValueError("Thematic data not loaded")
        geom_union = buffer_pos(
            self.aligner.get_thematic_union(), GRB_MAX_REFERENCE_BUFFER
        )
        collection = get_collection_grb_fiscal_parcels(
            year=self.year,
            geometry=geom_union,
            partition=self.part,
            crs=self.aligner.CRS,
        )
        self.input = dict(collection)
        self.aligner.logger.feedback_info(f"Adpf downloaded for year: {self.year}")
        return super().load_data()


class GRBSpecificDateParcelLoader(GeoJsonLoader):
    def __init__(self, date, aligner, partition=1000):
        logging.warning(
            "Loader for GRB parcel-situation on specific date (experimental); Use it with care!!!"
        )
        try:
            date = datetime.strptime(date, DATE_FORMAT).date()
            if date.year >= datetime.now().year:
                raise ValueError(
                    "The GRBSpecificDateParcelLoader can only be used for dates prior to the current year."
                )
        except Exception:
            raise ValueError(
                "No valid date, please provide a date in the format: " + DATE_FORMAT
            )
        super().__init__(_input=None, id_property=GRB_PARCEL_ID)
        self.aligner = aligner
        self.date = date
        self.part = partition
        self.data_dict_source["source"] = "Adp"
        self.data_dict_source[VERSION_DATE] = date.strftime(DATE_FORMAT)
        self.versiondate_info = {"name": GRB_VERSION_DATE, "format": datetime_format_TZ}

    def load_data(self):
        if not self.aligner.dict_thematic:
            raise ValueError("Thematic data not loaded")
        geom_union = buffer_pos(
            self.aligner.get_thematic_union(), GRB_MAX_REFERENCE_BUFFER
        )
        collection = get_collection_grb_parcels_by_date(
            date=self.date,
            geometry=geom_union,
            partition=self.part,
            crs=self.aligner.CRS,
        )
        self.input = dict(collection)
        self.aligner.logger.feedback_info(
            f"Parcels downloaded for specific date: {self.date.strftime(DATE_FORMAT)}"
        )
        return super().load_data()
