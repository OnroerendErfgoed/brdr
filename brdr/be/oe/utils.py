import logging

from shapely import box

from brdr.be.oe.enums import OEType
from brdr.constants import DOWNLOAD_LIMIT, DEFAULT_CRS
from brdr.geometry_utils import to_crs, from_crs
from brdr.utils import get_collection_by_partition

log = logging.getLogger(__name__)


def get_collection_oe_objects(
    oetype=OEType.AO,
    objectids=None,
    bbox=None,
    limit=DOWNLOAD_LIMIT,
    partition=-1,
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
        id_property = "uri"
    elif oetype == OEType.EO:
        typename = "lu:lu_wet_erfgobj_pub"
        id_property = "uri"
    else:
        logging.warning(
            "Undefined OE-type: " + str(oetype) + ": Empty collection returned"
        )
        return {}, None
    crs = to_crs(crs)
    theme_url = (
        "https://www.mercator.vlaanderen.be/raadpleegdienstenmercatorpubliek/wfs?"
    )
    params = {
        "SERVICE": "WFS",
        "VERSION": "2.0.0",
        "REQUEST": "GetFeature",
        "TYPENAMES": typename,
        "SRSNAME": from_crs(crs),
        "outputFormat": "application/json",
        "limit": limit,
    }
    if objectids is not None:
        params["CQL_FILTER"] = (
            id_property
            + " IN "
            + "("
            + ",".join([f"'{str(o)}'" for o in objectids])
            + ")"
        )
    bbox_polygon = None
    if bbox is not None:
        bbox_polygon = box(*tuple(o for o in bbox))
    collection = get_collection_by_partition(
        url=theme_url,
        params=params,
        geometry=bbox_polygon,
        partition=partition,
        crs=crs,
    )
    return (
        collection,
        id_property,
    )
