import json
import logging
from copy import deepcopy
from datetime import date
from datetime import datetime

import numpy as np
from shapely import intersects, Polygon
from shapely.geometry import shape
from shapely.geometry.base import BaseGeometry

from brdr.aligner import Aligner
from brdr.constants import (
    DEFAULT_CRS,
    LAST_VERSION_DATE,
    DATE_FORMAT,
    VERSION_DATE,
    FORMULA_FIELD_NAME,
)
from brdr.constants import DOWNLOAD_LIMIT
from brdr.constants import GRB_BUILDING_ID
from brdr.constants import GRB_FEATURE_URL
from brdr.constants import GRB_FISCAL_PARCELS_URL
from brdr.constants import GRB_KNW_ID
from brdr.constants import GRB_PARCEL_ID
from brdr.constants import GRB_VERSION_DATE
from brdr.constants import MAX_REFERENCE_BUFFER
from brdr.enums import GRBType
from brdr.geometry_utils import buffer_pos, safe_intersection
from brdr.geometry_utils import create_donut
from brdr.geometry_utils import features_by_geometric_operation
from brdr.geometry_utils import get_bbox
from brdr.loader import GeoJsonLoader, DictLoader
from brdr.logger import Logger
from brdr.utils import geojson_to_dicts
from brdr.utils import get_collection
from brdr.utils import get_collection_by_partition
from brdr.utils import get_series_geojson_dict

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


def get_geoms_affected_by_grb_change(
    aligner,
    grb_type=GRBType.ADP,
    date_start=date.today(),
    date_end=date.today(),
    one_by_one=False,
    border_distance=0,
):
    """
    Get a dictionary of thematic geometries that are affected bij GRB-changes in a
    specific timespan

    Args:
        aligner: Aligner instance
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
    dict_thematic = aligner.dict_thematic
    # if aligner.multi_as_single_modus:
    #     dict_thematic = merge_dict(dict_thematic)
    crs = aligner.CRS
    affected_dict: dict[str, BaseGeometry] = {}
    unchanged_dict: dict[str, BaseGeometry] = {}
    if border_distance > 0:
        for key in dict_thematic.keys():
            dict_thematic[key] = create_donut(dict_thematic[key], border_distance)
    if one_by_one:
        for key in dict_thematic:
            geom = dict_thematic[key]
            if is_grb_changed(geom, grb_type, date_start, date_end):
                affected_dict[key] = geom
            else:
                unchanged_dict[key] = geom
        return affected_dict, unchanged_dict
    else:
        # Temporal filter on VERDATUM
        geometry = aligner.get_thematic_union()
        coll_changed_grb, name_reference_id = get_collection_grb_actual(
            geometry,
            grb_type=grb_type,
            partition=1000,
            date_start=date_start,
            date_end=date_end,
            crs=crs,
        )
        dict_changed_grb, dict_changed_grb_properties = geojson_to_dicts(
            coll_changed_grb, name_reference_id
        )

        if len(dict_changed_grb) == 0:
            logging.info("No detected changes")
            return affected_dict, dict_thematic  # empty affected dict
        logging.info("Changed parcels in timespan: " + str(len(dict_changed_grb)))
        thematic_intersections = features_by_geometric_operation(
            list(dict_thematic.values()),
            list(dict_thematic.keys()),
            list(dict_changed_grb.values()),
            predicate="intersects",
        )
        logging.info("Number of filtered features: " + str(len(thematic_intersections)))
        for key, geom in dict_thematic.items():
            if key in thematic_intersections:
                affected_dict[key] = geom
            else:
                unchanged_dict[key] = geom
    return affected_dict, unchanged_dict


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
        GRB_FEATURE_URL + "/" + grb_type.upper() + "/items?"
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
        + grb_type.upper()
        + "/items?limit="
        + str(limit)
        + "&crs="
        + crs
    )
    if grb_type == GRBType.ADP:
        name_reference_id = GRB_PARCEL_ID
    elif grb_type == "gbg":
        name_reference_id = GRB_BUILDING_ID
    elif grb_type == GRBType.KNW:
        name_reference_id = GRB_KNW_ID
    else:
        logging.warning(
            f"type not implemented: {str(grb_type)} -->No reference-data loaded"
        )
        return {}, None

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
    formula_field=FORMULA_FIELD_NAME,
    max_distance_for_actualisation=2,
    feedback=None,
):
    """
    Function to update a thematic featurecollection to the most actual version of GRB.
    Important to notice that the featurecollection needs a 'formula' for the base-alignment.
    """
    logger = Logger(feedback)
    # Load featurecollection into a shapely_dict:
    dict_thematic = {}
    dict_thematic_props = {}

    last_version_date = datetime.now().date()
    for feature in featurecollection["features"]:
        id_theme = feature["properties"][id_theme_fieldname]
        try:
            geom = shape(feature["geometry"])
        except Exception:
            geom = Polygon()
        logger.feedback_debug("id theme: " + id_theme)
        logger.feedback_debug("geometry (wkt): " + geom.wkt)
        dict_thematic[id_theme] = geom
        try:
            dict_thematic_props[id_theme] = {
                FORMULA_FIELD_NAME: json.loads(feature["properties"][formula_field])
            }
            logger.feedback_debug("formula: " + str(dict_thematic_props[id_theme]))
        except Exception:
            raise Exception("Formula -attribute-field (json) cannot be loaded")
        try:
            logger.feedback_debug(str(dict_thematic_props[id_theme]))
            if (
                LAST_VERSION_DATE in dict_thematic_props[id_theme][FORMULA_FIELD_NAME]
                and dict_thematic_props[id_theme][FORMULA_FIELD_NAME][LAST_VERSION_DATE]
                is not None
                and dict_thematic_props[id_theme][FORMULA_FIELD_NAME][LAST_VERSION_DATE]
                != ""
            ):
                str_lvd = dict_thematic_props[id_theme][FORMULA_FIELD_NAME][
                    LAST_VERSION_DATE
                ]
                lvd = datetime.strptime(str_lvd, DATE_FORMAT).date()
                if lvd < last_version_date:
                    last_version_date = lvd
        except:
            raise Exception(f"Problem with {LAST_VERSION_DATE}")

    datetime_start = last_version_date
    datetime_end = datetime.now().date()
    base_aligner_result = Aligner(feedback=feedback)
    base_aligner_result.load_thematic_data(DictLoader(dict_thematic))
    base_aligner_result.name_thematic_id = id_theme_fieldname

    dict_affected, dict_unchanged = get_geoms_affected_by_grb_change(
        base_aligner_result,
        grb_type=GRBType.ADP,
        date_start=datetime_start,
        date_end=datetime_end,
        one_by_one=False,
    )
    logger.feedback_info(
        "Number of possible affected OE-thematic during timespan: "
        + str(len(dict_affected))
    )
    if len(dict_affected) == 0:
        logger.feedback_info(
            "No change detected in referencelayer during timespan. Script is finished"
        )
        return {}
    logger.feedback_debug(str(datetime_start))
    logger.feedback_debug(str(formula_field))

    # Initiate a Aligner to reference thematic features to the actual borders
    actual_aligner = Aligner(feedback=feedback)
    actual_aligner.load_thematic_data(
        DictLoader(data_dict=dict_affected, data_dict_properties=dict_thematic_props)
    )
    actual_aligner.load_reference_data(
        GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=actual_aligner)
    )

    actual_aligner.relevant_distances = (
        np.arange(0, max_distance_for_actualisation * 100, 10, dtype=int) / 100
    )
    dict_evaluated, prop_dictionary = actual_aligner.compare(
        threshold_area=5, threshold_percentage=1, dict_unchanged=dict_unchanged
    )

    return get_series_geojson_dict(
        dict_evaluated,
        crs=actual_aligner.CRS,
        id_field=actual_aligner.name_thematic_id,
        series_prop_dict=prop_dictionary,
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
        geom_union = buffer_pos(self.aligner.get_thematic_union(), MAX_REFERENCE_BUFFER)
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
        geom_union = buffer_pos(self.aligner.get_thematic_union(), MAX_REFERENCE_BUFFER)
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
        logging.warning("experimental loader; use with care!!!")
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
        geom_union = buffer_pos(self.aligner.get_thematic_union(), MAX_REFERENCE_BUFFER)
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
