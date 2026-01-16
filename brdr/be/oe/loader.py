import logging
from datetime import datetime

from brdr.be.oe.enums import OEType
from brdr.be.oe.utils import get_collection_oe_objects
from brdr.constants import DOWNLOAD_LIMIT, DEFAULT_CRS, DATE_FORMAT, VERSION_DATE
from brdr.geometry_utils import to_crs
from brdr.loader import GeoJsonLoader

log = logging.getLogger(__name__)


class OnroerendErfgoedLoader(GeoJsonLoader):
    def __init__(
        self,
        objectids=None,
        oetype=OEType.AO,
        bbox=None,
        limit=DOWNLOAD_LIMIT,
        partition=-1,
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
        self.crs = to_crs(crs)
        self.data_dict_source["source"] = "Onroerend Erfgoed"
        self.data_dict_source["source_url"] = "https://inventaris.onroerenderfgoed.be/"

    def load_data(self):
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
        logging.debug(f"OnroerendErfgoed-objects downloaded")
        return super().load_data()
