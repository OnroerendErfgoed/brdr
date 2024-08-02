import logging
from datetime import date, datetime

from shapely.geometry.base import BaseGeometry

from brdr.constants import (
    DOWNLOAD_LIMIT,
    DEFAULT_CRS,
    MAX_REFERENCE_BUFFER,
    GRB_FEATURE_URL,
)
from brdr.enums import GRBType
import numpy as np
from shapely import STRtree, intersects, prepare
from shapely.geometry import shape

from brdr.geometry_utils import create_dictionary_from_url
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
    # TODO add something that only detects changes in the border of objects, so big objects with internal parcel-updates are not affected?
    last_version_date = get_last_version_date(geometry, grb_type=grb_type)
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
):
    """
    Get thematic geometries that are affected bij GRB-changes in a specific timespan

    Args:
        dict_thematic: dictionary of thematic geometries
        grb_type: Type of GRB: parcels, buildings,...
        date_start: start-date to check changes in GRB
        date_end: end-date to check changes in GRB
        one_by_one: parameter to choose the methodology to check changes:
            * True: Every thematic geometry is checked individually (loop)
            * False: All GRB-parcels intersecting the thematic dictionary is checked at once
    Returns:
        dictionary/list of affected geometries

    """

    limit = DOWNLOAD_LIMIT
    crs = DEFAULT_CRS
    if one_by_one:
        affected_dict = {}
        for key in dict_thematic:
            geom = dict_thematic[key]
            if is_grb_changed(geom, grb_type, date_start, date_end):
                affected_dict[key] = geom
        return affected_dict
    else:
        # TODO: possible optimisation (if necessary) by using seperate calls for all features with smaller bbox (partitioning)

        # Temporal filter on VERDATUM
        # Example:https://geo.api.vlaanderen.be/GRB/ogc/features/v1/collections/ADP/items?limit=10&filter=VERSDATUM%20%3C%20%272023-07-01%27%20AND%20VERSDATUM%20%3E%20%272022-07-01%27&filter-lang=cql-text
        # Other filter capabilities on OGC feature API: https://portal.ogc.org/files/96288

        dict_changed_parcels, name_reference_id = get_reference_data_dict_grb_actual(
            dict_thematic=dict_thematic,
            grb_type=grb_type,
            crs=crs,
            partition=0,
            date_start=date_start,
            date_end=date_end,
        )
        affected_dict: dict[str, BaseGeometry] = {}

        # array_features = []
        # if "features" not in collection:
        #     logging.info("No detected changes")
        #     return affected_dict  # empty dict
        # for feature in collection["features"]:
        #     geom = shape(feature["geometry"]).buffer(0)
        #     array_features.append(geom)
        if len(dict_changed_parcels) == 0:
            logging.info("No detected changes")
            return affected_dict  # empty dict
        logging.info("Changed parcels in timespan: " + str(len(dict_changed_parcels)))
        # TODO extract this STRTREE-intersects to seperate 'geometry_utils'-function (also used in aligner)
        thematic_tree = STRtree(list(dict_thematic.values()))
        thematic_items = np.array(list(dict_thematic.keys()))
        arr_indices = thematic_tree.query(
            list(dict_changed_parcels.values()), predicate="intersects"
        )
        thematic_intersections = thematic_items.take(arr_indices[1])
        thematic_intersections = list(set(thematic_intersections))
        logging.info("Number of affected features: " + str(len(thematic_intersections)))
        for key in thematic_intersections:
            affected_dict[key] = dict_thematic[key]
    return affected_dict


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
    # TODO: performance optimisation? would it be possible and be faster to directly do a spatial filter (intersect), possible by feature API:
    # example? :https://geo.api.vlaanderen.be/GRB/ogc/features/collections/ADP/items?limit=100&crs=EPSG:31370&%20filter-crs=http://www.opengis.net/def/crs/EPSG/0/31370&filter-lang=cql-text&filter=INTERSECTS(SHAPE,Polygon%20((174210.85811228273087181%20172120.52806779573438689,%20174269.98297141725197434%20172110.67392460664268583,%20174280.9320194051542785%20172026.36625510000158101,%20174215.23773147788597271%20172027.46115989878308028,%20174210.85811228273087181%20172120.52806779573438689)))
    # https://portal.ogc.org/files/96288

    date_format = "%Y-%m-%d"
    limit = limit
    crs = crs
    bbox = str(geometry.bounds).strip("()")
    actual_url = (
        GRB_FEATURE_URL + "/" + grb_type.upper() + "/items?"
        "limit=" + str(limit) + "&crs=" + crs + "&bbox-crs=" + crs + "&bbox=" + bbox
    )
    update_dates = []
    collection = get_collection(actual_url, limit)
    if "features" not in collection:
        return None

    # This commented code seems slower then the uncommented version
    # dict_geoms = {}
    # dict_versdatum = {}
    # for c in collection["features"]:
    #     dict_geoms[c["properties"]["CAPAKEY"]] = shape(c["geometry"])
    #     dict_versdatum[c["properties"]["CAPAKEY"]] = c["properties"]["VERSDATUM"]
    # thematic_tree = STRtree(list(dict_geoms.values()))
    # thematic_items = np.array(list(dict_geoms.keys()))
    # arr_indices = thematic_tree.query([geometry], predicate="intersects")
    # thematic_intersections = thematic_items.take(arr_indices[1])
    # thematic_intersections = list(set(thematic_intersections))
    # logging.info("Number of intersected features: " + str(len(thematic_intersections)))
    # for key in thematic_intersections:
    #     versiondate = datetime.strptime(dict_versdatum[key], date_format).date()
    #     update_dates.append(versiondate)

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
    Fetches reference data (administrative plots, buildings, or artwork) from the GRB
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
    if not (date_start is None and date_end is None):
        date_format = "%Y-%m-%d"
        url_grb = (
            url_grb
            + "filter=VERSDATUM>"
            + date_start.strftime(date_format)
            + " AND VERSDATUM<"
            + date_end.strftime(date_format)
            + "&filter-lang=cql-text"
        )
    # url_grb = url_grb
    # + "&limit="
    # + str(limit)
    # + "&crs="
    # + crs
    # + "&bbox-crs="
    # + crs
    # + "&bbox="
    # + bbox
    limit = DOWNLOAD_LIMIT
    collection = {}
    bounds_array = []

    # Get the bounds of the thematic_data to get the necessary GRB-data
    for key in dict_thematic:
        # buffer theme geometry with x m (default 10)
        buffer_value = relevant_distance + MAX_REFERENCE_BUFFER

        # params: buffer_value, url,dict_thematic
        collection = create_dictionary_from_url(
            bounds_array,
            buffer_value,
            collection,
            crs,
            dict_thematic[key],
            key,
            limit,
            name_reference_id,
            partition,
            url_grb,
        )

    return collection, name_reference_id
