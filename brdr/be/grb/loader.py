import logging
from datetime import datetime

from brdr.be.grb.constants import GRB_MAX_REFERENCE_BUFFER, GRB_SUPPORTED_CRS
from brdr.be.grb.constants import GRB_PARCEL_ID
from brdr.be.grb.constants import GRB_VERSION_DATE
from brdr.be.grb.enums import GRBType
from brdr.be.grb.utils import (
    get_collection_grb_actual,
    get_collection_grb_fiscal_parcels,
    get_collection_grb_parcels_by_date,
)
from brdr.constants import (
    DATE_FORMAT,
    VERSION_DATE,
)
from brdr.geometry_utils import buffer_pos, to_crs
from brdr.loader import GeoJsonLoader

log = logging.getLogger(__name__)


datetime_format_TZ = "%Y-%m-%dT%H:%M:%SZ"

def check_crs(aligner):
    if not aligner.CRS in (to_crs(element) for element in GRB_SUPPORTED_CRS):
        raise ValueError(
            f"This GRB Loader only supports alignment in CRS '{GRB_SUPPORTED_CRS}' while CRS '{aligner.CRS}' is used"
        )


class GRBActualLoader(GeoJsonLoader):
    def __init__(self, grb_type: GRBType, aligner, partition: int = 1000):
        super().__init__()
        self.aligner = aligner
        check_crs(self.aligner)
        self.grb_type = grb_type
        self.part = partition
        self.data_dict_source["source"] = grb_type.value
        self.versiondate_info = {"name": GRB_VERSION_DATE, "format": DATE_FORMAT}



    def load_data(self):
        if not self.aligner.dict_thematic:
            raise ValueError("Thematic data not loaded")
        geom_union = buffer_pos(
            self.aligner.thematic_data.union, GRB_MAX_REFERENCE_BUFFER
        )
        collection, id_property = get_collection_grb_actual(
            grb_type=self.grb_type,
            geometry=geom_union,
            partition=self.part,
            crs=self.aligner.CRS,
        )
        self.id_property = id_property
        self.input = dict(collection)
        self.data_dict_source[VERSION_DATE] = datetime.now().strftime(DATE_FORMAT)
        self.aligner.logger.feedback_info(f"GRB downloaded: {self.grb_type}")
        return super().load_data()


class GRBFiscalParcelLoader(GeoJsonLoader):
    def __init__(self, year: str, aligner, partition=1000):
        super().__init__(_input=None, id_property=GRB_PARCEL_ID)
        self.aligner = aligner
        check_crs(self.aligner)
        self.year = year
        self.part = partition
        self.data_dict_source["source"] = "Adpf"
        self.data_dict_source[VERSION_DATE] = datetime(int(year), 1, 1).strftime(
            DATE_FORMAT
        )
        self.versiondate_info = {"name": GRB_VERSION_DATE, "format": datetime_format_TZ}

    def load_data(self):
        if not self.aligner.dict_thematic:
            raise ValueError("Thematic data not loaded")
        geom_union = buffer_pos(
            self.aligner.thematic_data.union, GRB_MAX_REFERENCE_BUFFER
        )
        collection = get_collection_grb_fiscal_parcels(
            year=self.year,
            geometry=geom_union,
            partition=self.part,
            crs=self.aligner.CRS,
        )
        self.input = dict(collection)
        self.aligner.logger.feedback_info(f"Adpf downloaded for year: {self.year}")
        return super().load_data()


class GRBSpecificDateParcelLoader(GeoJsonLoader):
    def __init__(self, date, aligner, partition=1000):
        logging.warning(
            "Loader for GRB parcel-situation on specific date (experimental); Use it with care!!!"
        )
        try:
            date = datetime.strptime(date, DATE_FORMAT).date()
            if date.year >= datetime.now().year:
                raise ValueError(
                    "The GRBSpecificDateParcelLoader can only be used for dates prior to the current year."
                )
        except Exception:
            raise ValueError(
                "No valid date, please provide a date in the format: " + DATE_FORMAT
            )
        super().__init__(_input=None, id_property=GRB_PARCEL_ID)
        self.aligner = aligner
        check_crs(self.aligner)
        self.date = date
        self.part = partition
        self.data_dict_source["source"] = "Adp"
        self.data_dict_source[VERSION_DATE] = date.strftime(DATE_FORMAT)
        self.versiondate_info = {"name": GRB_VERSION_DATE, "format": datetime_format_TZ}

    def load_data(self):
        if not self.aligner.dict_thematic:
            raise ValueError("Thematic data not loaded")
        geom_union = buffer_pos(
            self.aligner.thematic_data.union, GRB_MAX_REFERENCE_BUFFER
        )
        collection = get_collection_grb_parcels_by_date(
            date=self.date,
            geometry=geom_union,
            partition=self.part,
            crs=self.aligner.CRS,
        )
        self.input = dict(collection)
        self.aligner.logger.feedback_info(
            f"Parcels downloaded for specific date: {self.date.strftime(DATE_FORMAT)}"
        )
        return super().load_data()
