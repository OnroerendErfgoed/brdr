from datetime import datetime

from brdr.constants import (
    DATE_FORMAT,
    VERSION_DATE,
    DOWNLOAD_LIMIT,
)
from brdr.geometry_utils import buffer_pos, to_crs, from_crs
from brdr.loader import GeoJsonLoader
from brdr.nl.constants import (
    BRK_VERSION_DATE,
    BRK_FEATURE_URL,
    BRK_CRS,
    BRK_GENERIC_ID,
    BRK_MAX_REFERENCE_BUFFER,
)
from brdr.nl.enums import BRKType
from brdr.utils import get_collection_by_partition


def get_collection_brk(
    geometry,
    brk_type=BRKType.perceel,
    partition=1000,
    limit=DOWNLOAD_LIMIT,
    crs=BRK_CRS,
):
    crs = to_crs(crs)
    url = BRK_FEATURE_URL + "/" + brk_type.name + "/items?"

    name_reference_id = BRK_GENERIC_ID
    params = {"limit": limit, "crs": from_crs(crs), "f": "json"}

    collection = get_collection_by_partition(
        url=url,
        params=params,
        geometry=geometry,
        partition=partition,
        crs=from_crs(crs),
    )
    return collection, name_reference_id


class BRKLoader(GeoJsonLoader):
    def __init__(self, brk_type: BRKType, aligner, partition: int = 1000):
        super().__init__()
        self.aligner = aligner
        self.brk_type = brk_type
        self.part = partition
        self.data_dict_source["source"] = brk_type.value
        self.data_dict_source["url"] = BRK_FEATURE_URL + "/" + brk_type.name
        self.versiondate_info = {"name": BRK_VERSION_DATE, "format": DATE_FORMAT}

    def load_data(self):
        if not self.aligner.thematic_data:
            raise ValueError("Thematic data not loaded")
        if self.aligner.crs != to_crs(BRK_CRS):
            raise ValueError(
                f"BRKLoader only supports alignment in CRS '{BRK_CRS}' while CRS '{self.aligner.crs}' is used"
            )
        geom_union = buffer_pos(
            self.aligner.thematic_data.union, BRK_MAX_REFERENCE_BUFFER
        )
        if geom_union is None or geom_union.is_empty:
            raise ValueError(
                "Reference could not be loaded. Please load thematic data first"
            )

        collection, id_property = get_collection_brk(
            brk_type=self.brk_type,
            geometry=geom_union,
            partition=self.part,
            crs=self.aligner.crs,
        )
        self.id_property = id_property
        self.input = dict(collection)
        self.data_dict_source[VERSION_DATE] = datetime.now().strftime(DATE_FORMAT)
        self.aligner.logger.feedback_info(f"BRK downloaded: {self.brk_type}")
        return super().load_data()
