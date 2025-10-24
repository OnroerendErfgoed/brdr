from datetime import datetime

from brdr.constants import (
    DATE_FORMAT,
    BRK_VERSION_DATE,
    VERSION_DATE,
    DOWNLOAD_LIMIT,
    BRK_FEATURE_URL,
    BRK_CRS,
    BRK_GENERIC_ID,
    BRK_MAX_REFERENCE_BUFFER,
)
from brdr.geometry_utils import buffer_pos
from brdr.loader import GeoJsonLoader
from brdr.nl.enums import BRKType
from brdr.utils import get_collection_by_partition


def get_collection_brk(
    geometry,
    brk_type=BRKType.perceel,
    partition=1000,
    limit=DOWNLOAD_LIMIT,
    crs=BRK_CRS
):
    if crs == BRK_CRS:
        crs = 'http://www.opengis.net/def/crs/0/28992'
    else:
        raise ValueError (f"CRS expected: {BRK_CRS}, got CRS {crs} instead")
    url = (
        BRK_FEATURE_URL
        + "/"
        + brk_type.name
        + "/items?limit="
        + str(limit)
        + "&crs="
        + crs
        + "&f=json"
    )

    name_reference_id = BRK_GENERIC_ID
    collection = get_collection_by_partition(
        url, geometry=geometry, partition=partition, limit=limit, crs=crs
    )
    return collection, name_reference_id


class BRKLoader(GeoJsonLoader):
    def __init__(self, brk_type: BRKType, aligner, partition: int = 1000):
        super().__init__()
        self.aligner = aligner
        self.brk_type = brk_type
        self.part = partition
        self.data_dict_source["source"] = brk_type.value
        self.versiondate_info = {"name": BRK_VERSION_DATE, "format": DATE_FORMAT}

    def load_data(self):
        if not self.aligner.dict_thematic:
            raise ValueError("Thematic data not loaded")
        #TODO CRS support for other?
        if self.aligner.CRS!= BRK_CRS:
            raise ValueError(f"BRKLoader only supports alignment in CRS '{BRK_CRS}' while CRS '{self.aligner.CRS}' is used")
        geom_union = buffer_pos(
            self.aligner.get_thematic_union(), BRK_MAX_REFERENCE_BUFFER
        )

        collection, id_property = get_collection_brk(
            brk_type=self.brk_type,
            geometry=geom_union,
            partition=self.part,
            crs=self.aligner.CRS,
        )
        self.id_property = id_property
        self.input = dict(collection)
        self.data_dict_source[VERSION_DATE] = datetime.now().strftime(DATE_FORMAT)
        self.aligner.logger.feedback_info(f"BRK downloaded: {self.brk_type}")
        return super().load_data()
