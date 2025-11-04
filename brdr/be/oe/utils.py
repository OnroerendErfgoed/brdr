import logging

import requests
from shapely import box
from shapely.geometry import shape

from brdr.be.oe.enums import OEType
from brdr.constants import DOWNLOAD_LIMIT, DEFAULT_CRS
from brdr.utils import get_collection_by_partition

log = logging.getLogger(__name__)

def get_oe_dict_by_ids(objectids, oetype=OEType.AO):
    """
    Fetches thematic data for a list of objectIDs from the Inventaris Onroerend Erfgoed
    API.

    This function retrieves information about designated heritage objects
    (erfgoedobjecten or aanduidingsobjecten) from the Flemish Agency for Heritage (
    Inventaris Onroerend Erfgoed) based on a list of their IDs.

    Args:
        objectids (list): A list of objectIDs of 'erfgoedobjecten' or
            'aanduidingsobjecten'.
        oetype (string): A string: 'aanduidingsobjecten' (default) or 'erfgoedobjecten'

    Returns:
        dict: A dictionary where keys are objectIDs (as strings) and values are
              GeoJSON geometry objects. If an erfgoedobject/aanduidingsobject is not
              found, a corresponding warning message will be logged, but it won't be\
              included in the returned dictionary.

    Raises:
        requests.exceptions.RequestException: If there is an error fetching data from
            the API.
    """
    logging.warning("deprecated method, use OnroerendErfgoedLoader instead")
    # TODO remove function
    dict_thematic = {}
    if oetype == OEType.AO:
        typename = "aanduidingsobjecten"
        # id_property = "aanduid_id"
    elif oetype == OEType.EO:
        typename = "erfgoedobjecten"
        # id_property = "erfgoed_id"
    else:
        logging.warning(
            "Undefined OE-type: " + str(oetype) + ": Empty collection returned"
        )
        return {}, None

    base_url = "https://inventaris.onroerenderfgoed.be/" + typename + "/"
    headers = {"Accept": "application/json"}
    for a in objectids:
        url = base_url + str(a)
        response = requests.get(url, headers=headers).json()
        if "id" in response.keys():
            key = str(response["id"])
            geom = shape(response["locatie"]["contour"])
            dict_thematic[key] = geom
        else:
            logging.warning("object id " + str(a) + " not available in " + oetype)
    return dict_thematic


def get_collection_oe_objects(
    oetype=OEType.AO,
    objectids=None,
    bbox=None,
    limit=DOWNLOAD_LIMIT,
    partition=1000,
    crs=DEFAULT_CRS,
):
    """
    Fetches GeoJSON data for designated heritage objects (aanduidingsobjecten) within
    a bounding box.

    This function retrieves information about aanduidingsobjecten from the Flemish
    Mercator public WFS service using a bounding box (bbox) as a filter. The bbox should
    be provided in the format "xmin,ymin,xmax,ymax" (EPSG:31370 projection).

    Args:
        bbox (str): A comma-separated string representing the bounding box in EPSG:31370
                   projection (e.g., "100000,500000,200000,600000").
        limit (int, optional): The maximum number of features to retrieve per request.
                              Defaults to 1000.

    Returns:
        dict: A dictionary containing the retrieved GeoJSON feature collection. This
              collection might be truncated if the total number of features exceeds
              the specified limit.
    """
    if oetype == OEType.AO:
        typename = "ps:ps_aandobj"
        id_property = "aanduid_id"
    elif oetype == OEType.EO:
        typename = "lu:lu_wet_erfgobj_pub"
        id_property = "erfgoed_id"
    else:
        logging.warning(
            "Undefined OE-type: " + str(oetype) + ": Empty collection returned"
        )
        return {}, None

    theme_url ="https://www.mercator.vlaanderen.be/raadpleegdienstenmercatorpubliek/wfs?"
    params = {
        "SERVICE": "WFS",
        "VERSION": "2.0.0",
        "REQUEST": "GetFeature",
        "TYPENAMES": typename,
        "SRSNAME": crs,
        "outputFormat": "application/json",
        "limit": limit,
    }
    if objectids is not None:
        params["CQL_FILTER"] = f"{id_property} IN ("+ ", ".join(str(o) for o in objectids)+ ")"
    bbox_polygon = None
    if bbox is not None:
        bbox_polygon = box(*tuple(o for o in bbox))
    collection = get_collection_by_partition(
        url=theme_url,params=params, geometry=bbox_polygon, partition=partition, crs=crs
    )
    return (
        collection,
        id_property,
    )
