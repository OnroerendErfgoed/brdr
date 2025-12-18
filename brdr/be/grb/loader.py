import logging
from datetime import datetime
from typing import Any, Dict, Optional

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

def check_crs(aligner: Any) -> None:
    """
    Validates if the Aligner's Coordinate Reference System is supported by GRB.

    Args:
        aligner: The aligner object containing the CRS to check.

    Raises:
        ValueError: If the CRS used by the aligner is not in the supported GRB CRS list.
    """
    if not aligner.crs in (to_crs(element) for element in GRB_SUPPORTED_CRS):
        raise ValueError(
            f"This GRB Loader only supports alignment in CRS '{GRB_SUPPORTED_CRS}' while CRS '{aligner.crs}' is used"
        )


class GRBActualLoader(GeoJsonLoader):
    """
    Loader for the most recent version of GRB (Grootschalig Referentie Bestand) datasets.

    This loader fetches current thematic data (like buildings or roads) based on
    the spatial extent of the provided aligner.
    """

    def __init__(self, grb_type: GRBType, aligner: Any, partition: int = 1000):
        """
        Initialize the GRB Actual Loader.

        Args:
            grb_type: The specific type of GRB data to load (e.g., ADP, Gebouw).
            aligner: The aligner object providing the spatial context and CRS.
            partition: Number of features per request to handle large datasets.
                Defaults to 1000.
        """
        super().__init__()
        self.aligner = aligner
        check_crs(self.aligner)
        self.grb_type = grb_type
        self.part = partition
        self.data_dict_source["source"] = grb_type.value
        self.versiondate_info = {"name": GRB_VERSION_DATE, "format": DATE_FORMAT}

    def load_data(self) -> Any:
        """
        Downloads and processes the actual GRB data.

        Calculates a buffer around the thematic union of the aligner and
        fetches the requested GRB features.

        Returns:
            The result of the parent GeoJsonLoader's load_data method.

        Raises:
            ValueError: If thematic data has not been loaded into the aligner.
        """
        if not self.aligner.thematic_data:
            raise ValueError("Thematic data not loaded")
        geom_union = buffer_pos(
            self.aligner.thematic_data.union, GRB_MAX_REFERENCE_BUFFER
        )
        collection, id_property = get_collection_grb_actual(
            grb_type=self.grb_type,
            geometry=geom_union,
            partition=self.part,
            crs=self.aligner.crs,
        )
        self.id_property = id_property
        self.input = dict(collection)
        self.data_dict_source[VERSION_DATE] = datetime.now().strftime(DATE_FORMAT)
        self.aligner.logger.feedback_info(f"GRB downloaded: {self.grb_type}")
        return super().load_data()


class GRBFiscalParcelLoader(GeoJsonLoader):
    """
    Loader for fiscal parcel data (Adpf) for a specific year.
    """

    def __init__(self, year: str, aligner: Any, partition: int = 1000):
        """
        Initialize the Fiscal Parcel Loader.

        Args:
            year: The fiscal year to retrieve parcels for.
            aligner: The aligner object providing the spatial context.
            partition: Number of features per request. Defaults to 1000.
        """
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

    def load_data(self) -> Any:
        """
        Downloads and processes fiscal parcel data for the specified year.

        Returns:
            The result of the parent GeoJsonLoader's load_data method.

        Raises:
            ValueError: If thematic data has not been loaded into the aligner.
        """
        if not self.aligner.thematic_data:
            raise ValueError("Thematic data not loaded")
        geom_union = buffer_pos(
            self.aligner.thematic_data.union, GRB_MAX_REFERENCE_BUFFER
        )
        collection = get_collection_grb_fiscal_parcels(
            year=self.year,
            geometry=geom_union,
            partition=self.part,
            crs=self.aligner.crs,
        )
        self.input = dict(collection)
        self.aligner.logger.feedback_info(f"Adpf downloaded for year: {self.year}")
        return super().load_data()


class GRBSpecificDateParcelLoader(GeoJsonLoader):
    """
    Experimental loader for GRB parcel situations on a specific historical date.

    Note:
        This loader is intended for historical research and should be used with care.
    """

    def __init__(self, date: str, aligner: Any, partition: int = 1000):
        """
        Initialize the Specific Date Parcel Loader.

        Args:
            date: The historical date in the format defined by `DATE_FORMAT`.
            aligner: The aligner object providing spatial context.
            partition: Number of features per request. Defaults to 1000.

        Raises:
            ValueError: If the date is invalid or refers to the current year/future.
        """
        logging.warning(
            "Loader for GRB parcel-situation on specific date (experimental); Use it with care!!!"
        )
        try:
            date_obj = datetime.strptime(date, DATE_FORMAT).date()
            if date_obj.year >= datetime.now().year:
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
        self.date = date_obj
        self.part = partition
        self.data_dict_source["source"] = "Adp"
        self.data_dict_source[VERSION_DATE] = date_obj.strftime(DATE_FORMAT)
        self.versiondate_info = {"name": GRB_VERSION_DATE, "format": datetime_format_TZ}

    def load_data(self) -> Any:
        """
        Downloads and processes parcel data for the specified historical date.

        Returns:
            The result of the parent GeoJsonLoader's load_data method.

        Raises:
            ValueError: If thematic data has not been loaded into the aligner.
        """
        if not self.aligner.thematic_data:
            raise ValueError("Thematic data not loaded")
        geom_union = buffer_pos(
            self.aligner.thematic_data.union, GRB_MAX_REFERENCE_BUFFER
        )
        collection = get_collection_grb_parcels_by_date(
            date=self.date,
            geometry=geom_union,
            partition=self.part,
            crs=self.aligner.crs,
        )
        self.input = dict(collection)
        self.aligner.logger.feedback_info(
            f"Parcels downloaded for specific date: {self.date.strftime(DATE_FORMAT)}"
        )
        return super().load_data()