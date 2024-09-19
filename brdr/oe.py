import logging
from datetime import datetime
from enum import Enum

import requests
from shapely import box
from shapely.geometry import shape

from brdr.constants import DOWNLOAD_LIMIT, DEFAULT_CRS, DATE_FORMAT, VERSION_DATE
from brdr.loader import GeoJsonLoader
from brdr.logger import LOGGER
from brdr.utils import get_collection_by_partition

log = logging.getLogger(__name__)


class OEType(str, Enum):
    """
    Different types of Onroerend Eefgoed-objects are available:

    * AO: aanduidingsobjecten
    * EO: erfgoedobjecten
    """

    AO = "aanduidingsobjecten"
    EO = "erfgoedobjecten"


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

    theme_url = (
        "https://www.mercator.vlaanderen.be/raadpleegdienstenmercatorpubliek/wfs?"
        "SERVICE=WFS&REQUEST=GetFeature&VERSION=2.0.0&"
        f"TYPENAMES={typename}&"
        f"SRSNAME={crs}"
        "&outputFormat=application/json"
    )
    if objectids is not None:
        filter = (
            f"&CQL_FILTER={id_property} IN ("
            + ", ".join(str(o) for o in objectids)
            + ")"
        )
        theme_url = theme_url + filter
    bbox_polygon = None
    if bbox is not None:
        bbox_polygon = box(*tuple(o for o in bbox))

    return (
        get_collection_by_partition(
            theme_url, geometry=bbox_polygon, partition=partition, limit=limit, crs=crs
        ),
        id_property,
    )


class OnroerendErfgoedLoader(GeoJsonLoader):
    def __init__(
        self,
        objectids=None,
        oetype=OEType.AO,
        bbox=None,
        limit=DOWNLOAD_LIMIT,
        partition=1000,
        crs=DEFAULT_CRS,
    ):
        if (objectids is None and bbox is None) or (
            objectids is not None and bbox is not None
        ):
            raise ValueError("Please provide a ID-filter OR a BBOX-filter, not both")
        super().__init__()
        self.objectids = objectids
        self.oetype = oetype
        self.bbox = bbox
        self.limit = limit
        self.part = partition
        self.crs = crs
        self.data_dict_source["source"] = "Onroerend Erfgoed"

    def load_data(self):

        # geom_union = buffer_pos(self.aligner.get_thematic_union(), MAX_REFERENCE_BUFFER)
        collection, id_property = get_collection_oe_objects(
            oetype=self.oetype,
            objectids=self.objectids,
            bbox=self.bbox,
            partition=self.part,
            limit=self.limit,
            crs=self.crs,
        )
        self.id_property = id_property
        self.input = dict(collection)
        self.data_dict_source[VERSION_DATE] = datetime.now().strftime(DATE_FORMAT)
        LOGGER.debug(f"OnroerendErfgoed-objects downloaded")
        return super().load_data()
