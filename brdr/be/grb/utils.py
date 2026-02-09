import logging
from copy import deepcopy
from datetime import date
from datetime import datetime
from typing import Dict, Any

from shapely import intersects
from shapely.geometry import shape
from shapely.geometry.base import BaseGeometry

from brdr.be.grb.constants import (
    GRB_FEATURE_URL,
    GRB_FISCAL_PARCELS_URL,
    GRB_MAX_REFERENCE_BUFFER,
    GRB_PARCEL_ID,
    GRB_VERSION_DATE,
    GRB_GENERIC_ID,
)
from brdr.be.grb.enums import GRBType
from brdr.constants import DOWNLOAD_LIMIT, DEFAULT_CRS, DATE_FORMAT
from brdr.geometry_utils import (
    buffer_pos,
    safe_intersection,
    safe_unary_union,
    to_crs,
    from_crs,
)
from brdr.geometry_utils import create_donut
from brdr.geometry_utils import features_by_geometric_operation
from brdr.geometry_utils import get_bbox
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


def get_affected_ids_by_grb_change(
    thematic_geometries: Dict[Any, BaseGeometry],
    *,
    grb_type=GRBType.ADP,
    date_start=date.today(),
    date_end=date.today(),
    one_by_one=False,
    border_distance=0,
    geometry_thematic_union=None,
    crs=DEFAULT_CRS,
) -> list:
    """
    Get a list of affected IDs by GRB-changes in a
    specific timespan

    Args:
        thematic_geometries: dictionary if thematicID & Geometry
        grb_type: Type of GRB: parcels, buildings,...
        date_start: start-date to check changes in GRB
        date_end: end-date to check changes in GRB
        one_by_one: parameter to choose the methodology to check changes:
            * True: Every thematic geometry is checked individually (loop)
            * False: All GRB-parcels intersecting the thematic dictionary are checked
                at once
            (Default: False)
        border_distance: Distance that can be used to only check the 'border' of  the
            geometry, so 'big' geometries with internal parcel-updates are not affected
            (Default:0, indicating that the full geometry is checked for GRB-changes)
    Returns:
        list of affected IDs

    """

    affected = []
    if one_by_one:
        for key in thematic_geometries:
            geom = thematic_geometries[key]
            if is_grb_changed(
                geom, grb_type, date_start, date_end, border_distance=border_distance
            ):
                affected.append(key)
        return affected
    else:
        # Temporal filter on VERSIONDATE
        if geometry_thematic_union is None:
            geometry_thematic_union = safe_unary_union(
                list(thematic_geometries.values())
            )
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
            return []
        logging.info(
            f"Changed parcels in timespan with border_distance {str(border_distance)}: {str(len(dict_changed_grb))}"
        )
        thematic_intersections = features_by_geometric_operation(
            list(thematic_geometries.values()),
            list(thematic_geometries.keys()),
            list(dict_changed_grb.values()),
            predicate="intersects",
        )
        logging.info("Number of filtered features: " + str(len(thematic_intersections)))
        for key, geom in thematic_geometries.items():
            if key in thematic_intersections:
                affected.append(key)
    return affected


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
    crs = to_crs(crs)
    bbox = get_bbox(geometry)
    actual_url = GRB_FEATURE_URL + "/" + grb_type.name + "/items?"
    params = {
        "crs": from_crs(crs),
        "bbox-crs": from_crs(crs),
        "bbox": bbox,
        "limit": limit,
    }
    update_dates = []
    collection = get_collection(url=actual_url, params=params)
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
    crs = to_crs(crs)
    url = GRB_FEATURE_URL + "/" + grb_type.name + "/items?"
    params = {"crs": from_crs(crs), "limit": limit}

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
        params["filter"] = versiondate_filter
        params["filter-lang"] = "cql-text"

    collection = get_collection_by_partition(
        url=url, params=params, geometry=geometry, partition=partition, crs=crs
    )
    # TODO possibility to only return features that actually hit the unioned geometry, so the reference-data is not overloaded with features that are not relevant? Or maybe in a former/later step?
    return collection, name_reference_id


def get_collection_grb_fiscal_parcels(
    geometry,
    year=str(datetime.now().year),
    partition=1000,
    limit=DOWNLOAD_LIMIT,
    crs=DEFAULT_CRS,
):
    crs = to_crs(crs)
    url = GRB_FISCAL_PARCELS_URL + "/Adpf" + str(year) + "/items?"
    params = {"crs": from_crs(crs), "limit": limit}
    return get_collection_by_partition(
        url=url, params=params, geometry=geometry, partition=partition, crs=crs
    )


def get_collection_grb_parcels_by_date(
    geometry,
    date,
    partition=1000,
    limit=DOWNLOAD_LIMIT,
    crs=DEFAULT_CRS,
):
    crs = to_crs(crs)
    collection_year_after = get_collection_grb_fiscal_parcels(
        year=str(date.year),
        geometry=geometry,
        partition=partition,
        limit=limit,
        crs=crs,
    )
    # Filter on specific date: delete all features > specific_date
    # This is an experimental loader and results are not guaranteed; unclear if we have to use "year-1 & year" OR if we have to use "year & year + 1"
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
        # add them

    collection_specific_date = deepcopy(collection_year_after_filtered)
    filtered_features = collection_year_after_filtered["features"]
    specific_date_features = filtered_features + kept_features
    logging.debug(len(specific_date_features))
    collection_specific_date["features"] = specific_date_features

    return collection_specific_date
