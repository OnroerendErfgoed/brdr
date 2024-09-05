import json
import logging
from datetime import date
from datetime import datetime

from shapely import intersects
from shapely.geometry import shape
from shapely.geometry.base import BaseGeometry

from brdr.constants import DEFAULT_CRS
from brdr.constants import DOWNLOAD_LIMIT
from brdr.constants import GRB_BUILDING_ID
from brdr.constants import GRB_FEATURE_URL
from brdr.constants import GRB_FISCAL_PARCELS_URL
from brdr.constants import GRB_KNW_ID
from brdr.constants import GRB_PARCEL_ID
from brdr.constants import GRB_VERSION_DATE
from brdr.constants import MAX_REFERENCE_BUFFER
from brdr.enums import Evaluation
from brdr.enums import GRBType
from brdr.geometry_utils import buffer_pos
from brdr.geometry_utils import create_donut
from brdr.geometry_utils import features_by_geometric_operation
from brdr.geometry_utils import get_bbox
from brdr.loader import GeoJsonLoader
from brdr.utils import dict_series_by_keys
from brdr.utils import geojson_to_dicts
from brdr.utils import get_collection
from brdr.utils import get_collection_by_partition

log = logging.getLogger(__name__)
date_format = "%Y-%m-%d"


def is_grb_changed(
    geometry,
    grb_type=GRBType.ADP,
    date_start=date.today(),
    date_end=date.today(),
    border_distance=0,
    crs=DEFAULT_CRS,
):
    """
       checks if a geometry is possibly affected by changes in the reference layer
       during a specified timespan
    Args:
        geometry: Geometry to check on GRB-changes
        grb_type:  Type of GRB (parcels, buildings,artwork,...) to check
        date_start: Start of timespan to check if GRB-changes has occurred
        date_end: End of timespan to check if GRB-changes has occurred
        border_distance: Distance that can be used to only check the 'border' of
            the geometry, so 'big' geometries with internal parcel-updates are not
             affected
        (Default:0, indicating that the full geometry is checked fot GRB-changes)
        crs: Coordinate reference system to use

    Returns: Boolean, indicating if GRB is changed underneath the geometry

    """

    last_version_date = get_last_version_date(
        geometry, grb_type=grb_type, crs=crs, border_distance=border_distance
    )
    if last_version_date is None:
        return None
    if date_start <= last_version_date <= date_end:
        return True
    return False


def get_geoms_affected_by_grb_change(
    aligner,
    grb_type=GRBType.ADP,
    date_start=date.today(),
    date_end=date.today(),
    one_by_one=False,
    border_distance=0,
):
    """
    Get a dictionary of thematic geometries that are affected bij GRB-changes in a
    specific timespan

    Args:
        aligner: Aligner instance
        grb_type: Type of GRB: parcels, buildings,...
        date_start: start-date to check changes in GRB
        date_end: end-date to check changes in GRB
        one_by_one: parameter to choose the methodology to check changes:
            * True: Every thematic geometry is checked individually (loop)
            * False: All GRB-parcels intersecting the thematic dictionary is checked
                at once
        border_distance: Distance that can be used to only check the 'border' of  the
            geometry, so 'big' geometries with internal parcel-updates are not affected
            (Default:0, indicating that the full geometry is checked fot GRB-changes)
    Returns:
        dictionary of affected geometries

    """
    dict_thematic = aligner.dict_thematic
    # if aligner.multi_as_single_modus:
    #     dict_thematic = merge_dict(dict_thematic)
    crs = aligner.CRS
    affected_dict: dict[str, BaseGeometry] = {}
    unchanged_dict: dict[str, BaseGeometry] = {}
    if border_distance > 0:
        for key in dict_thematic.keys():
            dict_thematic[key] = create_donut(dict_thematic[key], border_distance)
    if one_by_one:
        for key in dict_thematic:
            geom = dict_thematic[key]
            if is_grb_changed(geom, grb_type, date_start, date_end):
                affected_dict[key] = geom
            else:
                unchanged_dict[key] = geom
        return affected_dict, unchanged_dict
    else:
        # Temporal filter on VERDATUM
        geometry = aligner.get_thematic_union()
        coll_changed_grb, name_reference_id = get_collection_grb_actual(
            geometry,
            grb_type=grb_type,
            partition=1000,
            date_start=date_start,
            date_end=date_end,
            crs=crs,
        )
        dict_changed_grb, dict_changed_grb_properties = geojson_to_dicts(
            coll_changed_grb, name_reference_id
        )

        if len(dict_changed_grb) == 0:
            logging.info("No detected changes")
            return affected_dict, dict_thematic  # empty affected dict
        logging.info("Changed parcels in timespan: " + str(len(dict_changed_grb)))
        thematic_intersections = features_by_geometric_operation(
            list(dict_thematic.values()),
            list(dict_thematic.keys()),
            list(dict_changed_grb.values()),
            predicate="intersects",
        )
        logging.info("Number of filtered features: " + str(len(thematic_intersections)))
        for key, geom in dict_thematic.items():
            if key in thematic_intersections:
                affected_dict[key] = geom
            else:
                unchanged_dict[key] = geom
    return affected_dict, unchanged_dict


def get_last_version_date(
    geometry,
    grb_type=GRBType.ADP,
    crs=DEFAULT_CRS,
    limit=DOWNLOAD_LIMIT,
    border_distance=0,
):
    """
    Retrieves the date of the last version for a specific geographic area within
    GRB (parcels, buildings,...).

    This function queries the GRB-API to find the most recent version-date (=last
    update  of object) for reference data of the specified `grb_type` (e.g., ADP,
    GBG, KNW) within the boundary of the provided `geometry`.

    Args:
    geometry (BaseGeometry): A Shapely geometry representing the area of interest.
    grb_type (GRBType, optional): The type of GRB data to consider.
        Defaults to GRBType.ADP (administrative plots).

    Returns:
    str: The date of the last version for the specified GRB data type within the area,
         formatted as a string according to the GRB API response (usually YYYY-MM-DD).

    None: If no data is found for the given geometry and GRB type combination.
    """
    if border_distance > 0:
        geometry = create_donut(geometry, border_distance)
    bbox = get_bbox(geometry)
    actual_url = (
        GRB_FEATURE_URL + "/" + grb_type.upper() + "/items?"
        "limit=" + str(limit) + "&crs=" + crs + "&bbox-crs=" + crs + "&bbox=" + bbox
    )
    update_dates = []
    collection = get_collection(actual_url, limit)
    if "features" not in collection:
        return None
    for c in collection["features"]:
        if intersects(geometry, shape(c["geometry"])):
            versiondate = datetime.strptime(
                c["properties"][GRB_VERSION_DATE], date_format
            ).date()
            update_dates.append(versiondate)

    update_dates = sorted(update_dates, reverse=True)
    if len(update_dates) > 0:
        return update_dates[0]
    return None


def get_collection_grb_actual(
    geometry,
    grb_type=GRBType.ADP,
    partition=1000,
    limit=DOWNLOAD_LIMIT,
    crs=DEFAULT_CRS,
    date_start=None,
    date_end=None,
):

    url = (
        GRB_FEATURE_URL
        + "/"
        + grb_type.upper()
        + "/items?limit="
        + str(limit)
        + "&crs="
        + crs
    )
    if grb_type == GRBType.ADP:
        name_reference_id = GRB_PARCEL_ID
    elif grb_type == "gbg":
        name_reference_id = GRB_BUILDING_ID
    elif grb_type == GRBType.KNW:
        name_reference_id = GRB_KNW_ID
    else:
        logging.warning(
            f"type not implemented: {str(grb_type)} -->No reference-data loaded"
        )
        return {}, None

    versiondate_filter = ""
    versiondate_filter_start = ""
    versiondate_filter_end = ""
    if date_start is not None:
        versiondate_filter_start = (
            GRB_VERSION_DATE + ">" + date_start.strftime(date_format)
        )
        versiondate_filter = versiondate_filter_start
    if date_end is not None:
        versiondate_filter_end = GRB_VERSION_DATE + "<" + date_end.strftime(date_format)
        versiondate_filter = versiondate_filter_end
    if not (date_start is None and date_end is None):
        versiondate_filter = versiondate_filter_start + " AND " + versiondate_filter_end
    if versiondate_filter != "":
        url = url + "&filter=" + versiondate_filter + "&filter-lang=cql-text"

    collection = get_collection_by_partition(
        url, geometry=geometry, partition=partition, limit=limit, crs=crs
    )
    return collection, name_reference_id


def get_collection_grb_fiscal_parcels(
    geometry,
    year=str(datetime.now().year),
    partition=1000,
    limit=DOWNLOAD_LIMIT,
    crs=DEFAULT_CRS,
):
    url = (
        GRB_FISCAL_PARCELS_URL + "/Adpf" + year + "/items?"
        "limit=" + str(limit) + "&crs=" + crs
    )
    return get_collection_by_partition(
        url, geometry=geometry, partition=partition, limit=limit, crs=crs
    )

def evaluate(
    actual_aligner,
    dict_series,
    dict_predicted,
    thematic_dict_formula,
    threshold_area=5,
    threshold_percentage=1,
    dict_unchanged=None,
):
    """
    evaluate affected geometries and give attributes to evaluate and decide if new
    proposals can be used
    """
    if dict_unchanged is None:
        dict_unchanged = {}
    theme_ids = list(dict_series_by_keys(dict_series).keys())
    dict_evaluated_result = {}
    prop_dictionary = {}
    # Fill the dictionary-structure with empty values
    for dist in dict_series.keys():
        dict_evaluated_result[dist] = {}
        prop_dictionary[dist] = {}
        for theme_id in theme_ids:
            prop_dictionary[dist][theme_id] = {}
        for theme_id in dict_unchanged.keys():
            prop_dictionary[dist][theme_id] = {}

    dict_predicted_keys = dict_series_by_keys(dict_predicted)

    for theme_id, dist_dict in dict_predicted_keys.items():
        equality = False
        for dist in sorted(dist_dict.keys()):
            if equality:
                break
            geomresult = dist_dict[dist][theme_id]["result"]
            actual_formula = actual_aligner.get_formula(geomresult)
            prop_dictionary[dist][theme_id]["formula"] = json.dumps(actual_formula)
            base_formula = None
            if theme_id in thematic_dict_formula:
                base_formula = thematic_dict_formula[theme_id]
            equality, prop = check_equality(
                base_formula,
                actual_formula,
                threshold_area,
                threshold_percentage,
            )
            if equality:
                dict_evaluated_result[dist][theme_id] = dict_predicted[dist][theme_id]
                prop_dictionary[dist][theme_id]["evaluation"] = prop
                break

    evaluated_theme_ids = list(dict_series_by_keys(dict_evaluated_result).keys())
    # fill where no equality is found/ The biggest predicted distance is returned as
    # proposal
    for theme_id in theme_ids:
        if theme_id not in evaluated_theme_ids:
            if len(dict_predicted_keys[theme_id].keys()) == 0:
                result = dict_series[0][theme_id]
                dict_evaluated_result[0][theme_id] = result
                prop_dictionary[0][theme_id]["formula"] = json.dumps(
                    actual_aligner.get_formula(result["result"])
                )
                prop_dictionary[0][theme_id]["evaluation"] = Evaluation.NO_PREDICTION_5
                continue
            #Add all predicted features so they can be manually checked
            for  dist in dict_predicted_keys[theme_id].keys():
                predicted_resultset = dict_predicted[dist][theme_id]
                dict_evaluated_result[dist][theme_id] = predicted_resultset
                prop_dictionary[dist][theme_id]["formula"] = json.dumps(
                    actual_aligner.get_formula(predicted_resultset["result"])
                )
                prop_dictionary[dist][theme_id][
                    "evaluation"
                ] = Evaluation.TO_CHECK_4

    for theme_id, geom in dict_unchanged.items():
        result = {"result": geom}
        dict_evaluated_result[0][theme_id] = result
        prop_dictionary[0][theme_id]["evaluation"] = Evaluation.NO_CHANGE_6
        prop_dictionary[0][theme_id]["formula"] = json.dumps(
            actual_aligner.get_formula(result["result"])
        )
    return dict_evaluated_result, prop_dictionary

def check_equality(
    base_formula, actual_formula, threshold_area=5, threshold_percentage=1
):
    """
    function that checks if 2 formulas are equal (determined by business-logic)
    """
    # TODO: research naar aanduid_id 116448 (equality na 0.5m), 120194 (1m)
    # TODO: research and implementation of following ideas
    # TODO: refine equality comparison, make it more generic
    # TODO: Add control of OD to equality-comparison (see case aanduid_id 120288)
    # ideas:
    # * If result_diff smaller than 0.x --> automatic update
    # * big polygons: If 'outer ring' has same formula (do net check inner side) -->
    #   automatic update
    # ** outer ring can be calculated: 1) negative buffer 2) original - buffered

    if base_formula is None or actual_formula is None:
        return False, Evaluation.NO_PREDICTION_5
    od_alike = False
    if base_formula["reference_od"] is None and actual_formula["reference_od"] is None:
        od_alike = True
    elif base_formula["reference_od"] is None or actual_formula["reference_od"] is None:
        od_alike = False
    elif (
                    (
                        abs(
                            base_formula["reference_od"]["area"]
                            - actual_formula["reference_od"]["area"]
                        )
                        * 100
                        / base_formula["reference_od"]["area"]
                    )
                    < threshold_percentage
                ):
        od_alike = True

    if (
        base_formula["reference_features"].keys()
        == actual_formula["reference_features"].keys() and od_alike
    ):
        if base_formula["full"] and base_formula["full"]:
            return True, Evaluation.EQUALITY_FORMULA_GEOM_1

        equal_reference_features = True
        for key in base_formula["reference_features"].keys():
            if (
                (
                    base_formula["reference_features"][key]["full"]
                    == actual_formula["reference_features"][key]["full"]
                )
                or (
                    abs(
                        base_formula["reference_features"][key]["area"]
                        - actual_formula["reference_features"][key]["area"]
                    )
                    > threshold_area
                )
                or (
                    (
                        abs(
                            base_formula["reference_features"][key]["area"]
                            - actual_formula["reference_features"][key]["area"]
                        )
                        * 100
                        / base_formula["reference_features"][key]["area"]
                    )
                    > threshold_percentage
                )
            ):
                equal_reference_features = False
        if equal_reference_features:
            return True, Evaluation.EQUALITY_FORMULA_2
    if base_formula["full"] and base_formula["full"] and od_alike:
        return True, Evaluation.EQUALITY_GEOM_3
    return False, Evaluation.NO_PREDICTION_5


class GRBActualLoader(GeoJsonLoader):
    def __init__(self, grb_type: GRBType, aligner, partition: int = 1000):
        super().__init__()
        self.aligner = aligner
        self.grb_type = grb_type
        self.part = partition
        self.data_dict_source["source"] = grb_type.value

    def load_data(self):
        if not self.aligner.dict_thematic:
            raise ValueError("Thematic data not loaded")
        geom_union = buffer_pos(self.aligner.get_thematic_union(), MAX_REFERENCE_BUFFER)
        collection, id_property = get_collection_grb_actual(
            grb_type=self.grb_type,
            geometry=geom_union,
            partition=self.part,
            crs=self.aligner.CRS,
        )
        self.id_property = id_property
        self.input = dict(collection)
        self.data_dict_source["version_date"] = datetime.now().strftime(date_format)
        self.aligner.logger.feedback_info(f"GRB downloaded: {self.grb_type}")
        return super().load_data()


class GRBFiscalParcelLoader(GeoJsonLoader):
    def __init__(self, year: str, aligner, partition=1000):
        super().__init__(_input=None, id_property=GRB_PARCEL_ID)
        self.aligner = aligner
        self.year = year
        self.part = partition
        self.data_dict_source["source"] = "Adpf"
        self.data_dict_source["version_date"] = datetime(int(year), 1, 1).strftime(
            date_format
        )

    def load_data(self):
        if not self.aligner.dict_thematic:
            raise ValueError("Thematic data not loaded")
        geom_union = buffer_pos(self.aligner.get_thematic_union(), MAX_REFERENCE_BUFFER)
        collection = get_collection_grb_fiscal_parcels(
            year=self.year,
            geometry=geom_union,
            partition=self.part,
            crs=self.aligner.CRS,
        )
        self.input = dict(collection)
        self.aligner.logger.feedback_info(f"Adpf downloaded for year: {self.year}")
        return super().load_data()
