import logging
from datetime import date, datetime

import numpy as np
from shapely.geometry.base import BaseGeometry
from brdr.constants import (
    DOWNLOAD_LIMIT,
    DEFAULT_CRS,
    MAX_REFERENCE_BUFFER,
    GRB_FEATURE_URL,
    GRB_FISCAL_PARCELS_URL,
)
from brdr.enums import GRBType
from shapely import intersects
from shapely.geometry import shape

from brdr.geometry_utils import (
    create_dictionary_from_url,
    features_by_geometric_operation,
    create_donut, get_bbox,
)
from brdr.utils import get_collection, dict_series_by_keys

log = logging.getLogger(__name__)
date_format = "%Y-%m-%d"


def is_grb_changed(
    geometry,
    grb_type=GRBType.ADP,
    date_start=date.today(),
    date_end=date.today(),
    border_distance=0,
    crs=DEFAULT_CRS,
):
    """
       checks if a geometry is possibly affected by changes in the reference layer during a specified timespan
    Args:
        geometry: Geometry to check on GRB-changes
        grb_type:  Type of GRB (parcels, buildings,artwork,...) to check
        date_start: Start of timespan to check if GRB-changes has occurred
        date_end: End of timespan to check if GRB-changes has occurred
        border_distance: Distance that can be used to only check the 'border' of the geometry, so 'big' geometries with internal parcel-updates are not affected
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
    dict_thematic,
    grb_type=GRBType.ADP,
    date_start=date.today(),
    date_end=date.today(),
    one_by_one=False,
    border_distance=0,
    crs=DEFAULT_CRS,
):
    """
    Get a dictionary of thematic geometries that are affected bij GRB-changes in a specific timespan

    Args:
        dict_thematic: dictionary of thematic geometries
        grb_type: Type of GRB: parcels, buildings,...
        date_start: start-date to check changes in GRB
        date_end: end-date to check changes in GRB
        one_by_one: parameter to choose the methodology to check changes:
            * True: Every thematic geometry is checked individually (loop)
            * False: All GRB-parcels intersecting the thematic dictionary is checked at once
        border_distance: Distance that can be used to only check the 'border' of the geometry, so 'big' geometries with internal parcel-updates are not affected
        (Default:0, indicating that the full geometry is checked fot GRB-changes)
    Returns:
        dictionary of affected geometries

    """
    if border_distance > 0:
        for key in dict_thematic.keys():
            dict_thematic[key] = create_donut(dict_thematic[key], border_distance)
    if one_by_one:
        affected_dict = {}
        for key in dict_thematic:
            geom = dict_thematic[key]
            if is_grb_changed(geom, grb_type, date_start, date_end):
                affected_dict[key] = geom
        return affected_dict
    else:
        # Temporal filter on VERDATUM
        dict_changed_parcels, name_reference_id = get_reference_data_dict_grb_actual(
            dict_thematic=dict_thematic,
            grb_type=grb_type,
            crs=crs,
            partition=0,
            date_start=date_start,
            date_end=date_end,
        )
        affected_dict: dict[str, BaseGeometry] = {}

        if len(dict_changed_parcels) == 0:
            logging.info("No detected changes")
            return affected_dict  # empty dict
        logging.info("Changed parcels in timespan: " + str(len(dict_changed_parcels)))
        thematic_intersections = features_by_geometric_operation(
            list(dict_thematic.values()),
            list(dict_thematic.keys()),
            list(dict_changed_parcels.values()),
            predicate="intersects",
        )
        logging.info("Number of filtered features: " + str(len(thematic_intersections)))
        for key in thematic_intersections:
            affected_dict[key] = dict_thematic[key]
    return affected_dict


def get_last_version_date(
    geometry,
    grb_type=GRBType.ADP,
    crs=DEFAULT_CRS,
    limit=DOWNLOAD_LIMIT,
    border_distance=0,
):
    """
    Retrieves the date of the last version for a specific geographic area within  GRB (parcels, buildings,...)).

    This function queries the GRB-API to find the most recent version-date (=last update of object)
    for reference data of the specified `grb_type` (e.g., ADP, GBG, KNW) within the boundary of the provided `geometry`.

    Args:
    geometry (BaseGeometry): A Shapely geometry representing the area of interest.
    grb_type (GRBType, optional): The type of GRB data to consider. Defaults to GRBType.ADP (administrative plots).

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
                c["properties"]["VERSDATUM"], date_format
            ).date()
            update_dates.append(versiondate)

    update_dates = sorted(update_dates, reverse=True)
    if len(update_dates) > 0:
        return update_dates[0]
    return None


def get_reference_data_dict_grb_actual(
    dict_thematic,
    relevant_distance=1,
    grb_type=GRBType.ADP,
    crs=DEFAULT_CRS,
    partition=0,
    date_start=None,
    date_end=None,
):
    """
    Fetches reference data (administrative plots, buildings, or artwork) from the actual GRB
    API based on thematic data.

    This function retrieves reference data from the Grootschalig Referentie
    Bestand (GRB) depending on the specified `grb_type` (e.g., administrative
    plots (ADP), buildings (GBG), or artwork (KNW)).
    It uses the bounding boxes of the geometries in the loaded thematic data
    (`self.aligner.dict_thematic`) to filter the relevant reference data
    geographically.

    Args:
        grb_type (GRBType, optional): The type of reference data to retrieve.
            Defaults to GRBType.ADP (administrative plots).
        partition (int, optional): If greater than zero, partitions the bounding box
            of the thematic data into a grid before fetching reference data by
            partition. Defaults to 0 (no partitioning).

    Returns:
        tuple: A tuple containing two elements:

            - dict: A dictionary where keys are reference data identifiers
              (as defined by `name_reference_id`) and values are GeoJSON  geometry
              objects representing the reference data.
            - str: The name of the reference data identifier property
              (e.g., "CAPAKEY" for ADP).

    Raises:
        ValueError: If an unsupported `grb_type` is provided.
    """

    url_grb = GRB_FEATURE_URL + "/" + grb_type.upper() + "/items?"
    if grb_type == GRBType.ADP:
        name_reference_id = "CAPAKEY"
    elif grb_type == "gbg":
        name_reference_id = "OIDN"
    elif grb_type == GRBType.KNW:
        name_reference_id = "OIDN"
    else:
        logging.warning(
            f"type not implemented: {str(grb_type)} -->No reference-data loaded"
        )
        return

    versiondate_filter = ""
    if date_start is not None:
        versiondate_filter_start = "VERSDATUM>" + date_start.strftime(date_format)
        versiondate_filter = versiondate_filter_start
    if date_end is not None:
        versiondate_filter_end = "VERSDATUM<" + date_end.strftime(date_format)
        versiondate_filter = versiondate_filter_end
    if not (date_start is None and date_end is None):
        versiondate_filter = versiondate_filter_start + " AND " + versiondate_filter_end
    if versiondate_filter != "":
        url_grb = url_grb + "filter=" + versiondate_filter + "&filter-lang=cql-text"

    limit = DOWNLOAD_LIMIT
    dictionary = {}
    bounds_array = []

    # Get the bounds of the thematic_data to get the necessary GRB-data
    for key in dict_thematic:
        # buffer theme geometry with x m (default 10)
        buffer_value = relevant_distance + MAX_REFERENCE_BUFFER

        dictionary = create_dictionary_from_url(
            bounds_array,
            buffer_value,
            dictionary,
            crs,
            dict_thematic[key],
            key,
            limit,
            name_reference_id,
            partition,
            url_grb,
        )

    return dictionary, name_reference_id


def get_collection_grb_fiscal_parcels(
    year=str(datetime.now().year), bbox=None, limit=DOWNLOAD_LIMIT, crs=DEFAULT_CRS
):
    # Load the Base reference data
    url = (
        GRB_FISCAL_PARCELS_URL + "/Adpf" + year + "/items?"
        "limit=" + str(limit) + "&crs=" + crs
    )
    if bbox is not None:
        url = url + "&bbox-crs=" + crs + "&bbox=" + bbox
    return get_collection(url, limit)

def evaluate (actual_aligner,thematic_dict_formula,series=np.arange(0, 200, 10, dtype=int) / 100,threshold_area=5,threshold_percentage=1):
    """
    evaluate affected geometries and give attributes to evaluate and decide if new proposals can be used
    """
    theme_ids = actual_aligner.dict_thematic.keys()
    dict_evaluated_result = {}
    prop_dictionary = {}
    for dist in series:
        dict_evaluated_result[dist]= {}
        prop_dictionary[dist] = {}
        for theme_id in theme_ids:
            prop_dictionary[dist] [theme_id] = {}

    dict_series,dict_predicted,diffs_dict= actual_aligner.predictor(series)

    for dist,dict_predicted_dist in dict_predicted.items():
        for theme_id in dict_predicted_dist.keys():
            results = dict_predicted_dist[theme_id]
            actual_formula = actual_aligner.get_formula(results["result"])
            prop_dictionary[dist][theme_id]["formula"] = actual_formula
            equality, property = check_equality(actual_aligner.dict_thematic[theme_id],thematic_dict_formula[theme_id], results["result"],actual_formula,threshold_area,threshold_percentage)
            if equality:
                dict_evaluated_result[dist][theme_id] = dict_predicted[dist][theme_id]
                prop_dictionary[dist][theme_id]["evaluation"] = property
                break

    dict_predicted_keys = dict_series_by_keys(dict_predicted)
    evaluated_theme_ids = list(dict_series_by_keys(dict_evaluated_result).keys())
    #fill where no equality is found/ The smallest predicted distance is returned as proposal
    for theme_id in theme_ids:
        if theme_id not in evaluated_theme_ids:
            if len(dict_predicted_keys[theme_id].keys())==0:
                dict_evaluated_result[0][theme_id] = dict_series[0][theme_id]
                prop_dictionary[0][theme_id]["evaluation"] = "NO_PREDICTION"
                continue
            smallest_predicted_dist = list(dict_predicted_keys[theme_id].keys())[0]
            dict_evaluated_result[smallest_predicted_dist][theme_id] = dict_predicted[smallest_predicted_dist][theme_id]
            prop_dictionary[smallest_predicted_dist][theme_id]["evaluation"] = "TO_CHECK"
    #add parameters if (base_formula["reference_features"][key]["full"] == actual_formula["reference_features"][key]["full"]) and (abs(base_formula["reference_features"][key]['area'] - actual_formula["reference_features"][key]['area']) > 5) and (abs(base_formula["reference_features"][key]['area'] - actual_formula["reference_features"][key]['area']) / base_formula[key]['area']> 0.05):
    return dict_evaluated_result, prop_dictionary

def check_equality(base_geometry, base_formula, actual_geometry, actual_formula,threshold_area=5,threshold_percentage=1):
    """
    function that checks if 2 formulas are equal (determined by business-logic)
    """
    # TODO: research naar aanduid_id 116448 (equality na 0.5m), 120194 (1m)
    # TODO: research and implementation of following ideas
    # TODO: refine equality comparison, make it more generic
    # TODO: Make Enum for equality
    # ideas:
    # * If result_diff smaller than 0.x --> automatic update
    # * big polygons: If 'outer ring' has same formula (do net check inner side) --> automatic update
    # ** outer ring can be calculated: 1) negative buffer 2) original - buffered

    if base_formula["reference_features"].keys()==actual_formula["reference_features"].keys():
        if base_formula["full"] and base_formula["full"]:
            return True, "EQUAL_FORMULA_FULL"
        for key in base_formula["reference_features"].keys():
            if (base_formula["reference_features"][key]["full"] == actual_formula["reference_features"][key]["full"]) and (abs(base_formula["reference_features"][key]['area'] - actual_formula["reference_features"][key]['area']) <threshold_area) and ((abs(base_formula["reference_features"][key]['area'] - actual_formula["reference_features"][key]['area'])*100 / base_formula["reference_features"][key]['area'])<threshold_percentage):
                return True, "EQUAL_FORMULA"
    if base_formula["full"] and base_formula["full"]:
        return True, "EQUAL_FULL"
    return False, "NO_EQUALITY_DETECTED"