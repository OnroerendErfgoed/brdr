import logging
from datetime import datetime
from typing import Any

from brdr.be.grb.constants import (
    GRB_MAX_REFERENCE_BUFFER,
    GRB_SUPPORTED_CRS,
    GRB_FEATURE_URL,
    GRB_PARCEL_ID,
    GRB_VERSION_DATE,
    GRB_FISCAL_PARCELS_URL,
)
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
    Validate if the Aligner's Coordinate Reference System is supported by GRB.

    The GRB (Grootschalig Referentie Bestand) service has specific CRS
    requirements for spatial queries.

    Parameters
    ----------
    aligner : Any
        The aligner object containing the CRS property to validate.

    Raises
    ------
    ValueError
        If the CRS used by the aligner is not present in the
        `GRB_SUPPORTED_CRS` list.
    """
    if not aligner.crs in (to_crs(element) for element in GRB_SUPPORTED_CRS):
        raise ValueError(
            f"This GRB Loader only supports alignment in CRS '{GRB_SUPPORTED_CRS}' "
            f"while CRS '{aligner.crs}' is used"
        )


class GRBActualLoader(GeoJsonLoader):
    """
    Loader for the most recent version of GRB (Grootschalig Referentie Bestand).

    This loader fetches live thematic reference data (such as buildings,
    parcels, or road edges) based on the spatial extent of the thematic
    data currently held by the aligner.

    Parameters
    ----------
    grb_type : GRBType
        The specific type of GRB data to load (e.g., ADP, Gebouw).
    aligner : Any
        The aligner object providing the spatial context, logger, and CRS.
    partition : int, optional
        Number of features per request to handle large datasets,
        by default 1000.

    Attributes
    ----------
    aligner : Any
        Reference to the parent aligner.
    grb_type : GRBType
        The requested GRB layer type.
    part : int
        The partitioning size for downloads.
    versiondate_info : dict
        Information about the versioning metadata format.
    """

    def __init__(self, grb_type: GRBType, aligner: Any, partition: int = 1000):
        super().__init__()
        self.aligner = aligner
        check_crs(self.aligner)
        self.grb_type = grb_type
        self.part = partition
        self.data_dict_source["source"] = grb_type.value
        self.data_dict_source["source_url"] = GRB_FEATURE_URL + "/" + grb_type.name
        self.versiondate_info = {"name": GRB_VERSION_DATE, "format": DATE_FORMAT}

    def load_data(self) -> Any:
        """
        Download and process the actual GRB data from the web service.

        Calculates a search buffer around the union of thematic geometries
        to ensure all relevant reference features are captured.

        Returns
        -------
        Any
            The processed FeatureCollection from the parent GeoJsonLoader.

        Raises
        ------
        ValueError
            If thematic data has not been loaded into the aligner prior
            to calling this method.
        """
        if not self.aligner.thematic_data:
            raise ValueError(
                "Reference could not be loaded. Please load thematic data first"
            )
        geom_union = buffer_pos(
            self.aligner.thematic_data.union, GRB_MAX_REFERENCE_BUFFER
        )
        if geom_union is None or geom_union.is_empty:
            raise ValueError(
                "Reference could not be loaded. Please load thematic data first"
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
        collection = super().load_data()
        for feature_id, feature in collection.features.items():
            if hasattr(feature, "data_id") and feature.data_id:
                match self.grb_type:
                    case GRBType.ADP:
                        feature.data_uri = f"https://data.vlaanderen.be/id/geometry/capakey/{feature.data_id}"
                    case GRBType.KNW:
                        feature.data_uri = f"https://data.vlaanderen.be/id/grb/kunstwerken/{feature.data_id}"
                    case GRBType.GBG:
                        feature.data_uri = f"https://data.vlaanderen.be/id/grb/gebouwen/{feature.data_id}"
        return collection


class GRBFiscalParcelLoader(GeoJsonLoader):
    """
    Loader for fiscal parcel data (Adpf) for a specific year.

    This loader retrieves the cadastral situation as it was registered
    for fiscal purposes at the start of a specific year.

    Parameters
    ----------
    year : str
        The fiscal year to retrieve parcels for (e.g., "2023").
    aligner : Any
        The aligner object providing spatial context and CRS.
    partition : int, optional
        Number of features per request, by default 1000.
    """

    def __init__(self, year: str, aligner: Any, partition: int = 1000):
        super().__init__(_input=None, id_property=GRB_PARCEL_ID)
        self.aligner = aligner
        check_crs(self.aligner)
        self.year = year
        self.part = partition
        self.data_dict_source["source"] = "Adpf"
        self.data_dict_source["source_url"] = (
            GRB_FISCAL_PARCELS_URL + "/Adpf" + str(year)
        )
        self.data_dict_source[VERSION_DATE] = datetime(int(year), 1, 1).strftime(
            DATE_FORMAT
        )
        self.versiondate_info = {"name": GRB_VERSION_DATE, "format": datetime_format_TZ}

    def load_data(self) -> Any:
        """
        Download and process fiscal parcel data for the specified year.

        Returns
        -------
        Any
            The result of the parent GeoJsonLoader's load_data method.

        Raises
        ------
        ValueError
            If thematic data is missing from the aligner.
        """
        if not self.aligner.thematic_data:
            raise ValueError(
                "Reference could not be loaded. Please load thematic data first"
            )
        geom_union = buffer_pos(
            self.aligner.thematic_data.union, GRB_MAX_REFERENCE_BUFFER
        )
        if geom_union is None or geom_union.is_empty:
            raise ValueError(
                "Reference could not be loaded. Please load thematic data first"
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
    Loader for GRB parcel situations on a specific historical date.

    This loader allows for high-precision historical reconstruction of parcel
    boundaries as they existed on a specific day.

    .. warning::
       This loader is experimental and intended for historical research.
       Data availability depends on the historical archives of the GRB provider.

    Parameters
    ----------
    date : str
        The historical date in the format defined by `DATE_FORMAT`.
    aligner : Any
        The aligner object providing spatial context.
    partition : int, optional
        Number of features per request, by default 1000.

    Raises
    ------
    ValueError
        - If the date format is invalid.
        - If the date refers to the current or a future year.
    """

    def __init__(self, date: str, aligner: Any, partition: int = 1000):
        logging.warning(
            "Loader for GRB parcel-situation on specific date (experimental); "
            "Use it with care!!!"
        )
        try:
            date_obj = datetime.strptime(date, DATE_FORMAT).date()
            if date_obj.year >= datetime.now().year:
                raise ValueError(
                    "The GRBSpecificDateParcelLoader can only be used for dates "
                    "prior to the current year."
                )
        except Exception:
            raise ValueError(
                f"No valid date, please provide a date in the format: {DATE_FORMAT}"
            )

        super().__init__(_input=None, id_property=GRB_PARCEL_ID)
        self.aligner = aligner
        check_crs(self.aligner)
        self.date = date_obj
        self.part = partition
        self.data_dict_source["source"] = "Adp"
        # This is a temporary source_url, that mimics the experimental implementation in
        # `get_collections_grb_parcels_by_date().
        self.data_dict_source["source_url"] = (
            GRB_FISCAL_PARCELS_URL + "/Adpf" + str(date_obj.year)
        )
        self.data_dict_source[VERSION_DATE] = date_obj.strftime(DATE_FORMAT)
        self.versiondate_info = {"name": GRB_VERSION_DATE, "format": datetime_format_TZ}

    def load_data(self) -> Any:
        """
        Download and process parcel data for the specified historical date.

        Returns
        -------
        Any
            The result of the parent GeoJsonLoader's load_data method.

        Raises
        ------
        ValueError
            If thematic data has not been loaded into the aligner.
        """
        if not self.aligner.thematic_data:
            raise ValueError(
                "Reference could not be loaded. Please load thematic data first"
            )

        geom_union = buffer_pos(
            self.aligner.thematic_data.union, GRB_MAX_REFERENCE_BUFFER
        )
        if geom_union is None or geom_union.is_empty:
            raise ValueError(
                "Reference could not be loaded. Please load thematic data first"
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
