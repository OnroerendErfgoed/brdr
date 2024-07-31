import logging
from datetime import date, datetime
from brdr.constants import DOWNLOAD_LIMIT, DEFAULT_CRS
from brdr.enums import GRBType
import numpy as np
from shapely import STRtree
from shapely.geometry import shape
from brdr.utils import get_collection

log = logging.getLogger(__name__)


def is_grb_changed(
    geometry,
    grb_type=GRBType.ADP,
    date_start=date.today(),
    date_end=date.today(),
):
    """
    checks if a geometry is possibly affected by changes in the reference layer during a specified timespan
    """
    last_version_date = get_last_version_date(geometry, grb_type=grb_type)
    if last_version_date is None:
        return None
    if last_version_date > date_start and last_version_date < date_end:
        return True
    return False


def get_geoms_affected_by_grb_change(
    dict_thematic,
    grb_type=GRBType.ADP,
    date_start=date.today(),
    date_end=date.today(),
    one_by_one=True,
):
    """
        # SEARCH FOR CHANGED Parcels in specific timespan
    # =================================================
    one_by_one: True. elke geometrie wordt beoordeeld
    False: We gaan globaal kijken in referentielaag waar wijzigingen zitten, en dan kijken en dan zien welke geoms geaffecteerd zijn
    """
    limit = DOWNLOAD_LIMIT
    crs = DEFAULT_CRS
    # bbox= te bepalen door dict
    if one_by_one:
        affected_dict = {}
        for key in dict_thematic:
            geom = dict_thematic[key]
            if is_grb_changed(geom, grb_type, date_start, date_end):
                affected_dict[key] = geom
        return affected_dict
    else:
        # TODO: check if temporal filter on VERSDATE is possible, by datetime, or by regex?: https://portal.ogc.org/files/96288
        bbox = []
        version_date = (
            "2023-01-01"  # all changes between this date and NOW will be found
        )
        base_url = (
            "https://geo.api.vlaanderen.be/GRB/ogc/features/collections/ADP/items?"
        )
        adp_url = (
            base_url
            + "datetime="
            + version_date
            + "%2F9999-12-31T00:00:00Z"
            + "& limit="
            + str(limit)
            + "&crs="
            + crs
            + "&bbox-crs=EPSG:31370&bbox="
            + bbox
        )
        print(adp_url)
        collection = get_collection(adp_url, limit)
        array_features = []
        for feature in collection["features"]:
            geom = shape(feature["geometry"]).buffer(0)
            array_features.append(geom)
        if len(array_features) == 0:
            logging.info("No detected changes")
            return {}
        logging.info("Changed parcels in timespan: " + str(len(array_features)))
        thematic_tree = STRtree(list(dict_thematic.values()))
        thematic_items = np.array(list(dict_thematic.keys()))
        arr_indices = thematic_tree.query(array_features, predicate="intersects")
        thematic_intersections = thematic_items.take(arr_indices[1])
        thematic_intersections = list(set(thematic_intersections))
        logging.info("Number of affected features: " + str(len(thematic_intersections)))
        # TODO create a dict as result

    return thematic_intersections


def get_last_version_date(
    geometry, grb_type=GRBType.ADP, crs=DEFAULT_CRS, limit=DOWNLOAD_LIMIT
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

    limit = limit
    crs = crs
    # TODO: Change implementation by bbox by a real intersection?: https://portal.ogc.org/files/96288
    # *directly possibly by feature API
    # *or, by filtering the bbox features on intersection
    bbox = str(geometry.bounds).strip("()")
    grb_type = grb_type.upper()
    actual_url = (
        "https://geo.api.vlaanderen.be/GRB/ogc/features/collections/"
        + grb_type
        + "/items?"
        "limit=" + str(limit) + "&crs=" + crs + "&bbox-crs=" + crs + "&bbox=" + bbox
    )
    update_dates = []
    collection = get_collection(actual_url, limit)
    if "features" not in collection:
        return None
    for c in collection["features"]:
        date_format = "%Y-%m-%d"
        versiondate = datetime.strptime(
            c["properties"]["VERSDATUM"], date_format
        ).date()
        update_dates.append(versiondate)

    update_dates = sorted(update_dates, reverse=True)

    return update_dates[0]
