import json
import logging
import os
from math import pi

import numpy as np
import requests as requests
from shapely import GeometryCollection, remove_repeated_points
from shapely import Polygon
from shapely import STRtree
from shapely import buffer
from shapely import equals
from shapely import get_exterior_ring
from shapely import get_interior_ring
from shapely import get_num_interior_rings
from shapely import get_parts
from shapely import make_valid
from shapely import polygons
from shapely import to_geojson
from shapely import unary_union
from shapely.geometry import shape
from shapely.geometry.base import BaseGeometry
from shapely.prepared import prep

from brdr.constants import BUFFER_MULTIPLICATION_FACTOR
from brdr.constants import CORR_DISTANCE
from brdr.constants import DOWNLOAD_LIMIT
from brdr.constants import DEFAULT_CRS
from brdr.constants import MAX_REFERENCE_BUFFER
from brdr.constants import MITRE_LIMIT
from brdr.constants import QUAD_SEGMENTS
from brdr.constants import THRESHOLD_CIRCLE_RATIO
from brdr.constants import THRESHOLD_EXCLUSION_AREA
from brdr.constants import THRESHOLD_EXCLUSION_PERCENTAGE
from brdr.enums import GRBType
from brdr.enums import OpenbaarDomeinStrategy
from brdr.geometry_utils import buffer_neg
from brdr.geometry_utils import buffer_neg_pos
from brdr.geometry_utils import buffer_pos
from brdr.geometry_utils import grid_bounds
from brdr.geometry_utils import safe_difference
from brdr.geometry_utils import safe_intersection
from brdr.geometry_utils import safe_symmetric_difference
from brdr.geometry_utils import safe_union
from brdr.utils import export_geojson, diffs_from_dict_series, get_breakpoints_zerostreak, \
    filter_resulting_series_by_key, get_collection

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(message)s", datefmt="%d-%b-%y %H:%M:%S"
)


###################


class Aligner:
    """
    This class is used to compare the thematic data with the reference data.
    The reference data can be loaded in different ways, for example by using the GRB data.
    The thematic data can be loaded by using a geojson file.
    The class can be used to compare the thematic data with the reference data.
    """

    def __init__(
        self,
        *,
        feedback=None,
        relevant_distance=1,
        threshold_overlap_percentage=50,
        od_strategy=OpenbaarDomeinStrategy.SNAP_SINGLE_SIDE,
        crs=DEFAULT_CRS
    ):
        """
        Initializes the Aligner object

        Args:
            feedback (object, optional): Feedback object that can be added to show
                feedback in QGIS. Defaults to None.
            relevant_distance (int, optional): The relevant distance (in meters) for
                processing. Defaults to 1.
            od_strategy (int, optional): The strategy to determine how to handle
                information outside the reference polygons (Openbaar Domein)
                (default 1: SNAP_SINGLE_SIDE)
            threshold_overlap_percentage (int, optional): Threshold (%) to determine
                from which overlapping-percentage a reference-polygon has to be included
                when there aren't relevant intersections or relevant differences
                (default 50%)

        """
        self.feedback = feedback
        self.relevant_distance = relevant_distance
        self.od_strategy = od_strategy
        self.threshold_overlap_percentage = threshold_overlap_percentage

        # PROCESSING DEFAULTS
        # thematic
        thematic_input = None  # to save the initially loaded geojson
        # name of the identifier-field of the thematic data (id has to be unique)
        self.name_thematic_id = "theme_identifier"
        # dictionary to store all thematic geometries to handle
        self.dict_thematic = {}

        # reference
        self.reference_input = None  # to save the initially loaded geojson

        self.name_reference_id = "ref_identifier" # name of the identifier-field of the reference data (id has to be unique,f.e CAPAKEY for GRB-parcels)
        self.dict_reference = {}  # dictionary to store all reference geometries
        self.reference_union = None # to save a unioned geometry of all reference polygons; needed for calculation in most OD-strategies

        # output-dictionaries (when processing dict_thematic)
        self.dict_result = None  # dictionary to save resulting geometries
        self.dict_result_diff = None  # dictionary to save global resulting differences
        self.dict_result_diff_plus = None  # dictionary to save positive resulting differences
        self.dict_result_diff_min = None  # dictionary to save negative resulting differences
        self.dict_relevant_intersection = None  # dictionary to save relevant_intersections
        self.dict_relevant_difference = None  # dictionary to save relevant_differences

        # Coordinate reference system
        # thematic geometries and reference geometries are assumed to be in the same CRS
        # before loading into the Aligner. No CRS-transformation will be performed.
        # When loading data, CRS is expected to be a projected CRS with units in 'meter (m)'.
        # By default EPSG:31370 (Lambert72), alternative: EPSG:3812 (Lambert2008)
        self.CRS = crs

        self.feedback_info("Aligner initialized")

    def feedback_debug(self, text):
        if self.feedback is not None:
            # self.feedback.pushInfo(text)
            return
        logging.debug(text)

    def feedback_info(self, text):
        if self.feedback is not None:
            self.feedback.pushInfo(text)
            return
        logging.info(text)

    def feedback_warning(self, text):
        if self.feedback is not None:
            self.feedback.pushInfo(text)
            return
        logging.warning(text)

    def buffer_distance(self):
        return self.relevant_distance / 2

    def load_thematic_data_geojson(self, thematic_input, name_thematic_id):
        self.thematic_input = thematic_input
        self.name_thematic_id = name_thematic_id
        self._prepare_thematic_data()

    def load_thematic_data_file(self, path_to_file, name_thematic_id):
        with open(path_to_file, "r") as f:
            thematic_input = json.load(f)
        self.load_thematic_data_geojson(thematic_input, name_thematic_id)

    def load_thematic_data_dict(self, dict_theme):
        self.dict_thematic = dict_theme

    def load_thematic_data_url(self, url, name_thematic_id):
        thematic_input = requests.get(url).json()
        self.load_thematic_data_geojson(thematic_input, name_thematic_id)

    def load_reference_data_dict(self, dict_ref):
        self.dict_reference = dict_ref
        self._prepare_reference_data()

    def load_reference_data_geojson(self, reference_input, name_reference_id):
        self.reference_input = reference_input
        self.name_reference_id = name_reference_id
        self._prepare_reference_data()

    def load_reference_data_file(self, path_to_file, name_reference_id):
        with open(path_to_file, "r") as f:
            reference_input = json.load(f)
        self.load_reference_data_geojson(reference_input, name_reference_id)

    def load_reference_data_url(self, url, name_reference_id):
        reference_input = requests.get(url).json()
        self.load_reference_data_geojson(reference_input, name_reference_id)

    def load_reference_data_grb_actual(self, *, grb_type=GRBType.ADP, partition=0):
        self.dict_reference = {}
        self.reference_input = None
        dict_ref, name_reference_id = self.get_reference_data_dict_grb_actual(
            grb_type, partition
        )
        self.name_reference_id = name_reference_id
        self.feedback_info(f"GRB downloaded: {grb_type}")
        self.load_reference_data_dict(dict_ref)

    def get_reference_data_dict_grb_actual(self, grb_type=GRBType.ADP, partition=0):
        """
        Fetches reference data (administrative plots, buildings, or artwork) from the GRB API based on thematic data.

        This function retrieves reference data from the Grootschalig Referentie Bestand (GRB)
        depending on the specified `grb_type` (e.g., administrative plots (ADP), buildings (GBG), or artwork (KNW)).
        It uses the bounding boxes of the geometries in the loaded thematic data (`self.dict_thematic`)
        to filter the relevant reference data geographically.

        Args:
            grb_type (GRBType, optional): The type of reference data to retrieve. Defaults to GRBType.ADP (administrative plots).
            partition (int, optional): If greater than zero, partitions the bounding box of the thematic data
                                       into a grid before fetching reference data by partition. Defaults to 0 (no partitioning).

        Returns:
            tuple: A tuple containing two elements:
                - dict: A dictionary where keys are reference data identifiers (as defined by `name_reference_id`)
                        and values are GeoJSON geometry objects representing the reference data.
                - str: The name of the reference data identifier property (e.g., "CAPAKEY" for ADP).

        Raises:
            ValueError: If an unsupported `grb_type` is provided.
        """
        if grb_type == GRBType.ADP:
            url_grb = (
                "https://geo.api.vlaanderen.be/GRB/ogc/features/collections/ADP/items?"
            )
            name_reference_id = "CAPAKEY"
        elif grb_type == "gbg":
            url_grb = (
                "https://geo.api.vlaanderen.be/GRB/ogc/features/collections/GBG/items?"
            )
            name_reference_id = "OIDN"
        elif grb_type == GRBType.KNW:
            url_grb = (
                "https://geo.api.vlaanderen.be/GRB/ogc/features/collections/KNW/items?"
            )
            name_reference_id = "OIDN"
        else:
            self.feedback_info(
                f"type not implemented: {str(grb_type)} -->No reference-data loaded"
            )
            return

        crs = self.CRS
        limit = DOWNLOAD_LIMIT
        collection = {}
        bounds_array = []

        # Get the bounds of the thematic_data to get the necessary GRB-data
        for key in self.dict_thematic:
            # buffer them geometry with x m (default 10)
            buffer_value = self.relevant_distance + MAX_REFERENCE_BUFFER
            geom = buffer(
                self.dict_thematic[key],
                buffer_value,
                quad_segs=QUAD_SEGMENTS,
                join_style="mitre",
                mitre_limit=MITRE_LIMIT,
            )
            bounds_array.append(geom)
            if partition < 1:
                bbox = str(geom.bounds).strip("()")
                url_grb_bbox = (
                    url_grb
                    + "f=application%2Fgeo%2Bjson&limit="
                    + str(limit)
                    + "&crs="
                    + crs
                    + "&bbox-crs="
                    + crs
                    + "&bbox="
                    + bbox
                )
                self.feedback_debug(key + "-->" + str(url_grb_bbox))
                coll = self._get_dict_from_url(url_grb_bbox, name_reference_id, limit)
                collection.update(coll)
            if partition > 0:
                geom = unary_union(bounds_array)
                grid = self.partition(geom, partition)
                for g in grid:
                    bbox = str(g.bounds).strip("()")
                    url_grb_bbox = (
                        url_grb
                        + "f=application%2Fgeo%2Bjson&limit="
                        + str(limit)
                        + "&crs="
                        + crs
                        + "&bbox-crs="
                        + crs
                        + "&bbox="
                        + bbox
                    )
                    self.feedback_debug(key + "-->" + str(url_grb_bbox))
                    coll = self._get_dict_from_url(
                        url_grb_bbox, name_reference_id, limit
                    )
                    collection.update(coll)

        return collection, name_reference_id

    def process_geometry(
        self,
        geometry: BaseGeometry,
        relevant_distance=1,
        od_strategy=OpenbaarDomeinStrategy.SNAP_SINGLE_SIDE,
        treshold_overlap_percentage=50,
    ) -> tuple[BaseGeometry, ...]:
        """
        method to align a geometry to the reference layer

        Args:
            geometry (BaseGeometry): The input geometric object.
            relevant_distance
            od_strategy
            treshold_overlap_percentage (float): The buffer distance (positive or negative).

        Returns:
            tuple[BaseGeometry...] : A tuple containing the resulting geometries:

            *   result (BaseGeometry): The resulting output geometry
            *   result_diff (BaseGeometry): The resulting difference output geometry
            *   result_diff_plus (BaseGeometry): The resulting positive difference output
                geometry
            *   result_diff_min (BaseGeometry): The resulting negative difference output
                geometry
            *   relevant_intersection (BaseGeometry): The relevant_intersection
            *   relevant_difference (BaseGeometry): The relevant_difference
        Notes:
            -
        Example:
        """
        self.feedback_debug("process geometry")
        self.relevant_distance = relevant_distance
        self.od_strategy = od_strategy
        self.threshold_overlap_percentage = treshold_overlap_percentage
        # array with all relevant parts of a thematic geometry; initial empty Polygon
        preresult = [Polygon()]
        (
            geometry,
            preresult,
            relevant_intersection_array,
            relevant_diff_array,
        ) = self._calculate_intersection_between_geometry_and_od(geometry, preresult)
        # get a list of all ref_ids that are intersecting the thematic geometry
        ref_intersections = self.reference_items.take(
            self.reference_tree.query(geometry)
        ).tolist()
        for key_ref in ref_intersections:
            geom_reference = self.dict_reference[key_ref]
            geom_intersection = safe_intersection(geometry, geom_reference)
            if geom_intersection.is_empty or geom_intersection is None:
                continue
            self.feedback_debug("calculate intersection")
            (
                geom,
                relevant_intersection,
                relevant_diff,
            ) = self._calculate_geom_by_intersection_and_reference(
                geom_intersection, geom_reference, False
            )
            self.feedback_debug("intersection calculated")
            preresult = self.add_multi_polygons_from_geom_to_array(geom, preresult)
            relevant_intersection_array = self.add_multi_polygons_from_geom_to_array(
                relevant_intersection, relevant_intersection_array
            )
            relevant_diff_array = self.add_multi_polygons_from_geom_to_array(
                relevant_diff, relevant_diff_array
            )
        # UNION INTERMEDIATE LAYERS
        relevant_intersection = unary_union(relevant_intersection_array)
        if relevant_intersection is None or relevant_intersection.is_empty:
            relevant_intersection = Polygon()
        relevant_diff = unary_union(relevant_diff_array)
        if relevant_diff is None or relevant_diff.is_empty:
            relevant_diff = Polygon()
        # POSTPROCESSING
        (
            result,
            result_diff,
            result_diff_plus,
            result_diff_min,
        ) = self._postprocess_preresult(preresult, geometry)

        return (
            result,
            result_diff,
            result_diff_plus,
            result_diff_min,
            relevant_intersection,
            relevant_diff,
        )

    def process_dict_thematic(
        self,
        relevant_distance=1,
        od_strategy=OpenbaarDomeinStrategy.SNAP_SINGLE_SIDE,
        treshold_overlap_percentage=50,
    ):
        """
        Aligns a thematic dictionary of geometries to the reference layer based on
        specified parameters. - method to align a thematic dictionary to the reference
        layer

        Args:
            relevant_distance (float, optional): The relevant distance (in meters) for
                processing. Defaults to 1.
            od_strategy (int, optional): The strategy for overlap detection.
                Defaults to 1.
            treshold_overlap_percentage (float, optional): The threshold percentage for
                considering full overlap. Defaults to 50.

        Returns:
            tuple: A tuple containing dictionaries with processed data:
                - dict_result: Aligned thematic data for each thematic key.
                - dict_result_diff: global differences between thematic data and reference data.
                - dict_result_diff_plus: Positive differences.
                - dict_result_diff_min: Negative differences.
                - dict_relevant_intersection: relevant intersections.
                - dict_relevant_diff: relevant differences.

        """

        dict_result = {}
        dict_result_diff = {}
        dict_result_diff_plus = {}
        dict_result_diff_min = {}
        dict_relevant_intersection = {}
        dict_relevant_diff = {}
        for key in self.dict_thematic:
            self.feedback_info("thematic id to process: " + str(key))
            (
                result,
                result_diff,
                result_diff_plus,
                result_diff_min,
                relevant_intersection,
                relevant_diff,
            ) = self.process_geometry(
                self.dict_thematic[key],
                relevant_distance,
                od_strategy,
                treshold_overlap_percentage,
            )
            dict_result[key] = result
            dict_result_diff[key] = result_diff
            dict_result_diff_plus[key] = result_diff_plus
            dict_result_diff_min[key] = result_diff_min
            dict_relevant_intersection[key] = relevant_intersection
            dict_relevant_diff[key] = relevant_diff
        self.dict_result = dict_result
        self.dict_result_diff = dict_result_diff
        self.dict_result_diff_plus = dict_result_diff_plus
        self.dict_result_diff_min = dict_result_diff_min
        self.dict_relevant_intersection = dict_relevant_intersection
        self.dict_relevant_difference = dict_relevant_diff
        self.feedback_info("thematic dictionary processed")
        return (
            dict_result,
            dict_result_diff,
            dict_result_diff_plus,
            dict_result_diff_min,
            dict_relevant_intersection,
            dict_relevant_diff,
        )

    def predictor(
            self,
            relevant_distances=np.arange(0, 300, 10, dtype=int)/100,
            od_strategy=OpenbaarDomeinStrategy.SNAP_SINGLE_SIDE,
            treshold_overlap_percentage=50,
    ):
        """
        Predicts the 'most interesting' relevant distances for changes in thematic elements based on a distance series.

        This function analyzes a set of thematic geometries (`self.dict_thematic`) to identify potentially 
        interesting distances where changes occur. It performs the following steps:

        1. **Process Distance Series:**
            - Calculates a series of results for different distances specified by `relevant_distances`.
            - This calculation might involve functions like `self.process_series` (implementation details likely depend on your specific code).

        2. **Calculate Difference Metrics:**
            - Analyzes the results from the distance series to compute difference metrics 
              between thematic elements at each distance (using `diffs_from_dict_series`).

        3. **Identify Breakpoints and Zero-Streaks:**
            - For each thematic geometry, it identifies potential "breakpoints" where the difference metric changes sign (from positive to negative or vice versa).
            - It also identifies "zero-streaks" which are consecutive distances with a difference metric close to zero (potentially indicating minimal change).

        4. **Predict Interesting Distances:**
            - The function considers distances corresponding to breakpoints and zero-streaks as potentially interesting for further analysis.
            - These distances are stored in a dictionary (`dict_predicted`) with the thematic element key as the outer key.
            - Additionally, the corresponding results (tuples) from the distance series for those distances are included.

        5. **Filter Results:**
            - The function might further filter the predicted results for each thematic element based on the element key (using `filter_resulting_series_by_key`).

        Args:
            relevant_distances (np.ndarray, optional): A NumPy array of distances to be analyzed. Defaults to np.arange(0.1, 5.05, 0.1).
            od_strategy (OpenbaarDomeinStrategy, optional): A strategy for handling open data in the processing (implementation specific). Defaults to OpenbaarDomeinStrategy.SNAP_SINGLE_SIDE.
            treshold_overlap_percentage (int, optional): A percentage threshold for considering full overlap in the processing (implementation specific). Defaults to 50.

        Returns:
            dict: A dictionary containing predicted interesting distances for each thematic element.
                - Keys: Thematic element identifiers from `self.dict_thematic`.
                - Values: Dictionaries with the following structure for each thematic element:
                    - Keys: Distances identified as interesting (breakpoints or zero-streaks).
                    - Values: Tuples containing results (likely specific to your implementation) from the distance series for the corresponding distance.

        Logs:
            - Debug logs the thematic element key being processed.
        """
        dict_predicted = {}
        for key in self.dict_thematic.keys():
            dict_predicted[key]={}
        dict_series = self.process_series(relevant_distances=relevant_distances,od_strategy=od_strategy,treshold_overlap_percentage=treshold_overlap_percentage)
        diffs = diffs_from_dict_series(dict_series, self.dict_thematic)
        for key in diffs:
            if len(diffs[key]) == len(relevant_distances):
                lst_diffs = list(diffs[key].values())
                breakpoints, zero_streaks = get_breakpoints_zerostreak(relevant_distances, lst_diffs)
                logging.debug(str(key))
                for zs in zero_streaks:
                    dict_predicted[key][zs[0]] = dict_series[zs[0]]
                dict_predicted[key] = filter_resulting_series_by_key(dict_predicted[key],key)
        return dict_predicted, diffs

    def process_series(
        self,
        relevant_distances,
        od_strategy=OpenbaarDomeinStrategy.SNAP_SINGLE_SIDE,
        treshold_overlap_percentage=50,
    ):
        """
        Calculates the resulting dictionaries for thematic data based on a series of relevant
            distances.

        Args:
            relevant_distances (list): A list of relevant distances (in meters) to
                process.
            od_strategy (int, optional): The strategy for overlap detection.
                Defaults to 1.
            treshold_overlap_percentage (float, optional): The threshold percentage for
                considering full overlap. Defaults to 50.

        Returns:
            dict: A dictionary containing the resulting dictionaries for a series of relevant distances:

                {
                    'relevant_distance_1': (tuple of resulting dictionaries),
                    'relevant_distance_2': (tuple of resulting dictionaries),
                    ...
                }
        """
        self.feedback_debug("Process series" + str(relevant_distances))
        self.od_strategy = od_strategy
        self.threshold_overlap_percentage = treshold_overlap_percentage
        #self._prepare_thematic_data() #not necessary? Assumed that dict_thematic is already loaded
        dict_series = {}
        for s in relevant_distances:
            self.feedback_info(
                "Processing series - relevant_distance (m):"
                + str(s)
                + " with ODStrategy "
                + str(self.od_strategy)
            )
            dict_series[s] = self.process_dict_thematic(s, od_strategy)
        self.feedback_info("End of processing series: " + str(relevant_distances))
        return dict_series

    def get_formula(self, geometry, with_geom=False):
        """
        Calculates formula-related information based on the input geometry.

        Args:
            geometry (shapely.geometry object): The input geometry.
            with_geom (bool, optional): Whether to include geometry information in the
                output. Defaults to False.

        Returns:
            dict: A dictionary containing formula-related data:

                -   'full': True if the intersection is the same as the reference
                    geometry, else False.
                -   'area': Area of the intersection or reference geometry.
                -   'percentage': Percentage of intersection area relative to the
                    reference geometry.
                -   'geometry': GeoJSON representation of the intersection (if
                    with_geom is True).
        """
        dict_formula = {}
        ref_intersections = self.reference_items.take(
            self.reference_tree.query(geometry)
        ).tolist()
        for key_ref in ref_intersections:
            geom_reference = self.dict_reference[key_ref]
            geom_intersection = make_valid(safe_intersection(geometry, geom_reference))
            if geom_intersection.is_empty or geom_intersection is None:
                continue
            if equals(geom_intersection, geom_reference):
                full = True
                area = round(geom_reference.area, 2)
                perc = 100
                geom = to_geojson(geom_reference)
            else:
                perc = round(geom_intersection.area * 100 / geom_reference.area, 2)
                if perc < 0.01:
                    continue
                elif perc > 99.99:
                    full = True
                    area = round(geom_reference.area, 2)
                    perc = 100
                    geom = to_geojson(geom_reference)
                else:
                    full = False
                    area = round(geom_intersection.area, 2)
                    geom = to_geojson(geom_intersection)
            if not with_geom:
                geom = None

            dict_formula[key_ref] = {
                "full": full,
                "area": area,
                "percentage": perc,
                "geometry": geom,
            }

        self.feedback_debug(str(dict_formula))
        return dict_formula

    def get_last_version_date(self, geometry, grb_type=GRBType.ADP):
        """
        Retrieves the date of the last version for a specific geographic area within  GRB (parcels, buildings,...)).

        This function queries the GRB-API to find the most recent version-date (=last update of object)
        for reference data of the specified `grb_type` (e.g., ADP, GBG, KNW) within the boundary of the provided `geometry`.

        Args:
            geometry (BaseGeometry): A Shapely geometry representing the area of interest.
            grb_type (GRBType, optional): The type of GRB data to consider. Defaults to GRBType.ADP (administrative plots).

        Returns:
            str: The date of the last version for the specified GRB data type within the area,
                 formatted as a string according to the GRB API response (usually YYYY-MM-DD).

            None: If no data is found for the given geometry and GRB type combination.
        """
        limit = DOWNLOAD_LIMIT
        crs = self.CRS
        bbox = str(geometry.bounds).strip("()")
        if grb_type is None:
            grb_type = "adp"
        grb_type = grb_type.upper()
        actual_url = (
            "https://geo.api.vlaanderen.be/GRB/ogc/features/collections/"
            + grb_type
            + "/items?"
            "limit=" + str(limit) + "&crs=" + crs + "&bbox-crs=" + crs +"&bbox=" + bbox
        )
        update_dates = []
        collection = get_collection(actual_url, limit)
        if "features" not in collection:
            return None
        for c in collection["features"]:
            update_dates.append(c["properties"]["VERSDATUM"])
        update_dates = sorted(update_dates, reverse=True)
        return update_dates[0]

    def export_results(self, path, multi_to_single=True):
        """
            Exports analysis results as GeoJSON files.

            This function exports 6 GeoJSON files containing the analysis results to the specified `path`.

            Args:
                path (str): The path to the directory where the GeoJSON files will be saved.
                multi_to_single (bool, optional): If True (default), converts MultiPolygon geometries to single Polygons
                                                    in the exported GeoJSON files. This can be useful for visualization purposes.

            Details of exported files:
                - result.geojson: Contains the original thematic data from `self.dict_result`.
                - result_diff.geojson: Contains the difference between the original and predicted data from `self.dict_result_diff`.
                - result_diff_plus.geojson: Contains results for areas that are added (increased area).
                - result_diff_min.geojson: Contains results for areas that are removed (decreased area).
            """
        export_geojson(
            os.path.join(path, "result.geojson"),
            self.dict_result,
            self.CRS,
            self.name_thematic_id,
            multi_to_single=multi_to_single
        )
        export_geojson(
            os.path.join(path, "result_diff.geojson"),
            self.dict_result_diff,
            self.CRS,
            self.name_thematic_id,
            multi_to_single=multi_to_single,
        )
        export_geojson(
            os.path.join(path, "result_diff_plus.geojson"),
            self.dict_result_diff_plus,
            self.CRS,
            self.name_thematic_id,
            multi_to_single=multi_to_single,

        )
        export_geojson(
            os.path.join(path, "result_diff_min.geojson"),
            self.dict_result_diff_min,
            self.CRS,
            self.name_thematic_id,
            multi_to_single=multi_to_single,
        )
        export_geojson(
            os.path.join(path, "relevant_intersection.geojson"),
            self.dict_relevant_intersection,
            self.CRS,
            self.name_thematic_id,
            multi_to_single=multi_to_single,
        )
        export_geojson(
            os.path.join(path, "relevant_difference.geojson"),
            self.dict_relevant_difference,
            self.CRS,
            self.name_thematic_id,
            multi_to_single=multi_to_single,
        )

    def _prepare_reference_data(self):
        """
        Prepares reference data for spatial queries and analysis.

        This function processes the reference data provided in `self.reference_input` (assumed to be a GeoJSON feature collection).
        It performs the following tasks:

        1. **Iterates through features:**
            - Extracts the reference data identifier (`key`) from each feature's properties using `self.name_reference_id`.
            - Converts the feature's geometry to a Shapely geometry object using `shape`.
            - Validates the geometry using `make_valid` to ensure it's a well-formed GeoJSON object.
            - Stores the processed data in `self.dict_reference` with the identifier as the key and the validated geometry as the value.
            - Provides debug feedback for each processed feature (key and geometry).

        2. **Optimizes spatial queries:**
            - Creates a Spatial Relationship Tree (STRtree) using `STRtree` for efficient spatial queries against the reference data in `self.dict_reference`.
            - Converts the dictionary keys (reference identifiers) to a NumPy array for potential performance benefits in future operations.

        3. **Clears reference union:**
            - Sets `self.reference_union` to `None`. This variable stores the combined geometry of all reference data,
              and it's cleared here to indicate that it needs to be recalculated if requested later.

        Returns:
            None
        """
        if self.reference_input is not None:
            for f in self.reference_input["features"]:
                key = f["properties"][self.name_reference_id]
                key = str(key)
                geom = shape(f["geometry"])
                self.dict_reference[key] = make_valid(geom)
                self.feedback_debug(key + "-->" + str(geom))
        # create an SRTree for performance optimisation
        self.feedback_info("length of reference_dict: " + str(len(self.dict_reference)))
        self.reference_tree = STRtree(list(self.dict_reference.values()))
        self.reference_items = np.array(list(self.dict_reference.keys()))
        # clear the reference_union, so it will be recalculated on request when needed
        self.reference_union = None
        return

    def _get_dict_from_url(self, input_url, name_reference_id, limit):
        collection = get_collection(input_url, limit)
        dictionary = {}
        if "features" not in collection or len(collection["features"]) == 0:
            return dictionary
        for f in collection["features"]:
            key = str(f["properties"][name_reference_id])
            geom = shape(f["geometry"])
            if key not in collection:
                dictionary[key] = make_valid(geom)
            self.feedback_debug(key + "-->" + str(geom))
        return dictionary

    def _calculate_intersection_between_geometry_and_od(self, geometry, preresult):
        # Calculate the intersection between thematic and Openbaar Domein
        relevant_intersection_array = [Polygon()]
        relevant_difference_array = [Polygon()]
        geom_thematic_od = Polygon()

        if self.od_strategy == OpenbaarDomeinStrategy.EXCLUDE:
            # Completely exclude everything that is not on the reference layer
            self.feedback_debug("OD-strategy EXCLUDE")
            # Remove from the thematic layer all parts that are not on the reference layer
            # !!this strategy adapts the input-geometry!!
            geometry = safe_intersection(geometry, self._get_reference_union())
        elif self.od_strategy == OpenbaarDomeinStrategy.AS_IS:
            # All parts that are not covered by the reference layer are added to the
            #         resulting geometry AS IS
            self.feedback_debug("OD-strategy AS IS")
            # all OD-parts wil be added AS IS
            geom_thematic_od = safe_difference(geometry, self._get_reference_union())
        elif self.od_strategy == OpenbaarDomeinStrategy.SNAP_SINGLE_SIDE:
            # Everything that falls within the relevant distance over the
            #  plot boundary is snapped to the plot.
            #  Only the inner-reference-boundaries are used.
            #  The outer-reference-boundaries are not used.
            self.feedback_debug("OD-strategy SNAP_SINGLE_SIDE")
            # geom of OD
            geom_od = safe_difference(geometry, self._get_reference_union())
            # only the relevant parts of OD
            geom_od_neg_pos = buffer_neg_pos(geom_od, self.buffer_distance())
            # geom_thematic_od = safe_intersection(geom_od_neg_pos,geom_od)# resulting
            # thematic OD
            geom_od_neg_pos_buffered = buffer_pos(
                geom_od_neg_pos, self.buffer_distance()
            )  # include parts
            geom_thematic_od = safe_intersection(
                geom_od_neg_pos_buffered, geom_od
            )  # resulting thematic OD
        elif self.od_strategy == OpenbaarDomeinStrategy.SNAP_ALL_SIDE:
            # Everything that falls within the relevant distance over
            # the plot boundary is snapped to the plot.
            #  Inner-reference-boundaries and outer-reference-boundaries are used.
            self.feedback_debug("OD-strategy SNAP BOTH SIDED")
            (
                geom_thematic_od,
                relevant_difference_array,
                relevant_intersection_array,
            ) = self._od_snap_all_side(
                geometry,
                relevant_difference_array,
                relevant_intersection_array,
            )
        elif self.od_strategy == OpenbaarDomeinStrategy.SNAP_FULL_AREA_SINGLE_SIDE:
            # Strategy useful for bigger areas.
            # integrates the entire inner area of the input geometry,
            # so Openbaar Domein of the inner area is included in the result
            # Combines SNAP_SINGLE_SIDE with the inner area
            self.feedback_debug("OD-strategy Full-area-variant of OD-SNAP_SINGLE_SIDE")
            geom_thematic_od = self._od_full_area(geometry)
        elif self.od_strategy == OpenbaarDomeinStrategy.SNAP_FULL_AREA_ALL_SIDE:
            # Strategy useful for bigger areas.
            # integrates the entire inner area of the input geometry,
            # so Openbaar Domein of the inner area is included in the result
            # Combines SNAP_ALL_SIDE with the inner area
            self.feedback_debug("OD-strategy Full-area-variant of OD-SNAP_ALL_SIDE")
            # first part is a copy of OD_ALL_SIDE
            (
                geom_thematic_od,
                relevant_difference_array,
                relevant_intersection_array,
            ) = self._od_snap_all_side(
                geometry,
                relevant_difference_array,
                relevant_intersection_array,
            )
            # This part is a copy of SNAP_FULL_AREA_SINGLE_SIDE
            geom_theme_od_min_clipped_plus_buffered_clipped = self._od_full_area(
                geometry
            )
            # UNION the calculation of  OD-SNAP_ALL_SIDE with FULL AREA of OD-SNAP_FULL_AREA_SINGLE_SIDE
            geom_thematic_od = safe_union(
                geom_theme_od_min_clipped_plus_buffered_clipped, geom_thematic_od
            )

        elif self.od_strategy == OpenbaarDomeinStrategy.SNAP_SINGLE_SIDE_VARIANT_1:
            # OD-strategy SNAP_SINGLE_SIDE - variant 1:
            # Everything that falls within the relevant distance over the
            #  plot boundary is snapped to the plot.
            #  Only the inner-reference-boundaries are used.
            #  The outer-reference-boundaries are not used.
            self.feedback_debug("OD-strategy SNAP_SINGLE_SIDE - variant 1")
            # geom of OD
            geom_od = safe_difference(geometry, self._get_reference_union())
            # only the relevant parts of OD
            geom_od_neg_pos = buffer_neg_pos(geom_od, self.buffer_distance())
            # resulting thematic OD
            geom_thematic_od = safe_intersection(geom_od_neg_pos, geom_od)
        elif self.od_strategy == OpenbaarDomeinStrategy.SNAP_SINGLE_SIDE_VARIANT_2:
            # OD-strategy SNAP_SINGLE_SIDE - variant 2:
            # Everything that falls within the relevant distance over the
            #  plot boundary is snapped to the plot.
            #  Only the inner-reference-boundaries are used.
            #  The outer-reference-boundaries are not used.
            self.feedback_debug("OD-strategy SNAP_SINGLE_SIDE - variant 2")
            # TODO: Drop?
            pass

        # ADD THEMATIC_OD
        preresult = self.add_multi_polygons_from_geom_to_array(
            geom_thematic_od, preresult
        )
        return (
            geometry,
            preresult,
            relevant_intersection_array,
            relevant_difference_array,
        )

    def _od_full_area(self, geometry):
        geom_theme_od = safe_difference(geometry, self._get_reference_union())
        geom_theme_min_buffered = buffer_neg(
            buffer_pos(
                buffer_neg(geometry, self.relevant_distance),
                self.buffer_distance(),
            ),
            self.buffer_distance(),
        )
        geom_theme_od_clipped_min_buffered = safe_intersection(
            geom_theme_min_buffered, geom_theme_od
        )
        geom_theme_od_min_clipped_plus_buffered = buffer_pos(
            geom_theme_od_clipped_min_buffered, self.relevant_distance
        )
        geom_theme_od_min_clipped_plus_buffered_clipped = safe_intersection(
            geom_theme_od_min_clipped_plus_buffered, geom_theme_od
        )
        geom_thematic_od = geom_theme_od_min_clipped_plus_buffered_clipped
        return geom_thematic_od

    def _od_snap_all_side(
        self,
        geometry,
        relevant_difference_array,
        relevant_intersection_array,
    ):
        geom_thematic_buffered = make_valid(
            buffer_pos(geometry, BUFFER_MULTIPLICATION_FACTOR * self.relevant_distance)
        )
        clip_ref_thematic_buffered = safe_intersection(
            self._get_reference_union(), geom_thematic_buffered
        )
        geom_reference = safe_difference(
            geom_thematic_buffered, clip_ref_thematic_buffered
        )  # beide OD-stukken worden SNAPPED toegevoegd
        geom_thematic = geometry
        if geom_reference.is_empty or geom_reference is None:
            geom_thematic_od = geom_reference
        else:
            geom_intersection = safe_intersection(geom_reference, geom_thematic)
            (
                geom_thematic_od,
                geom_relevant_intersection,
                geom_relevant_diff,
            ) = self._calculate_geom_by_intersection_and_reference(
                geom_intersection, geom_reference, True
            )
            relevant_intersection_array = self.add_multi_polygons_from_geom_to_array(
                geom_relevant_intersection, relevant_intersection_array
            )
            relevant_difference_array = self.add_multi_polygons_from_geom_to_array(
                geom_relevant_diff, relevant_difference_array
            )
        return geom_thematic_od, relevant_difference_array, relevant_intersection_array

    def _calculate_geom_by_intersection_and_reference(
        self,
        geom_intersection: BaseGeometry,
        geom_reference: BaseGeometry,
        is_openbaar_domein,
    ):
        """
        Calculates the geometry based on intersection and reference geometries.

        Args:
            geom_intersection (BaseGeometry): The intersection geometry.
            geom_reference (BaseGeometry): The reference geometry.
            is_openbaar_domein (bool): A flag indicating whether it's a public domain
                (area not covered with reference polygon).

        Returns:
            tuple: A tuple containing the resulting geometries:

            *   geom: BaseGeometry or None: The resulting geometry or None if conditions
                are not met.
            *   geom_relevant_intersection: BaseGeometry or None: The relevant
                intersection.
            *   geom_relevant_difference: BaseGeometry or None: The relevant difference.

        Notes:
            -   If the reference geometry area is 0, the overlap is set to 100%.
            -   If the overlap is less than relevant_OVERLAP_PERCENTAGE or the
                intersection area is less than relevant_OVERLAP_AREA, None is returned.
            -   Otherwise, the relevant intersection and difference geometries are
                calculated.
            -   If both relevant intersection and difference are non-empty, the final
                geometry is obtained by applying safe intersection and buffering.
            -   If only relevant intersection is non-empty, the result is the reference
                geometry.
            -   If only relevant difference is non-empty, the result is None.
        """
        if geom_reference.area == 0:
            overlap = 100
        else:
            overlap = geom_intersection.area * 100 / geom_reference.area
        if (
            overlap < THRESHOLD_EXCLUSION_PERCENTAGE
            or geom_intersection.area < THRESHOLD_EXCLUSION_AREA
        ):
            return Polygon(), Polygon(), Polygon()
        geom_difference = safe_difference(geom_reference, geom_intersection)
        geom_relevant_intersection = buffer_neg(
            geom_intersection, self.buffer_distance()
        )
        geom_relevant_difference = buffer_neg(geom_difference, self.buffer_distance())
        if (
            not geom_relevant_intersection.is_empty
            and not geom_relevant_difference.is_empty
        ):
            # intersectie en difference relevant
            geom_x = safe_intersection(
                geom_reference,
                safe_difference(
                    geom_reference,
                    safe_intersection(
                        geom_difference,
                        buffer_neg_pos(geom_difference, self.buffer_distance()),
                    ),
                ),
            )
            geom = safe_intersection(
                geom_x,
                buffer_pos(
                    buffer_neg_pos(geom_x, self.buffer_distance()),
                    self.buffer_distance(),
                ),
            )
            #TODO BEGIN: experimental fix - check if it is ok in all cases?
            #when calculating for OD, we create a 'virtual parcel'. When calculating this virtual parcel, it is buffered to take outer boundaries into account.
            #This results in a side-effect that there are extra non-logical parts included in the result. The function below tries to exclude these non-logica parts.
            # see eo_id 206363 with relevant distance=0.2m and SNAP_ALL_SIDE
            if is_openbaar_domein:
                #geom = buffer_neg_pos(geom, self.buffer_distance())
                geom = self.get_relevant_polygons_from_geom (geom)
            #TODO END
        elif (
            not geom_relevant_intersection.is_empty
            and geom_relevant_difference.is_empty
        ):
            geom = geom_reference
        elif (
            geom_relevant_intersection.is_empty
            and not geom_relevant_difference.is_empty
        ):
            # TODO: check needed
            # if overlap > treshold_overlap_percentage and openbaar domein:
            #     geom = snap_geom_to_reference(
            #       geom_intersection, geom_reference, relevant_distance
            #   )
            # else:
            geom = geom_relevant_intersection  # (=empty geometry)
        else:
            if is_openbaar_domein:
                geom = geom_relevant_intersection  # (=empty geometry)
            # geom = snap_geom_to_reference (geom_intersection, geom_reference,
            # relevant_distance)
            elif self.threshold_overlap_percentage < 0:
                # if we take a value of -1, the original border will be used
                geom = geom_intersection
            elif overlap > self.threshold_overlap_percentage:
                geom = geom_reference
            else:
                geom = geom_relevant_intersection  # (=empty geometry)
        return geom, geom_relevant_intersection, geom_relevant_difference


    def get_relevant_polygons_from_geom(self, geom):
        """
        Get only the relevant parts (polygon) from a geometry.
        Points, Lines and Polygons smaller than relevant distance are excluded from the result
        """
        if geom.is_empty or geom is None:
            # If the input geometry is empty or None, do nothing.
            return geom
        else:
            geom = make_valid(unary_union(geom))
            # Create a GeometryCollection from the input geometry.
            geometry_collection = GeometryCollection(geom)
            array=[]
            for g in geometry_collection.geoms:
                # Ensure each sub-geometry is valid.
                g = make_valid(g)
                if str(g.geom_type) in ["Polygon", "MultiPolygon"]:
                    relevant_geom = buffer_neg(g,self.buffer_distance())
                    if relevant_geom != None and not relevant_geom.is_empty:
                        array.append(g)
        return make_valid(unary_union(array))

    @staticmethod
    def _add_geom_to_dict(dictionary, geom, id_theme):
        dictionary[id_theme] = geom
        return

    # def _snap_geom_to_reference(self, geom_input, geom_reference, relevant_distance):
    # Deze functie werkt niet correct met Shapely. Deze vermijdt dat polygonen
    # collapsen als alles bijeen wordt gesnapt, wat we in sommige gevallen
    # effectief willen.
    # return snap(geom_input, geom_reference, relevant_distance)

    def _prepare_thematic_data(self):
        """
        Load geometries of a GeoJSON and stores them in a dictionary.

        This method processes the thematic data from the input GeoJSON file. It
        iterates through each feature, extracts the relevant properties, converts the
        geometry to a valid shape, and stores it in a dictionary.

        Returns:
            None.
        """
        # THEMATIC PREPARATION
        for f in self.thematic_input["features"]:
            key = f["properties"][self.name_thematic_id]
            key = str(key)
            geom = shape(f["geometry"])
            self.dict_thematic[key] = make_valid(geom)
            self.feedback_debug(key + "-->" + str(geom))
        return

    def _get_reference_union(self):
        if self.reference_union is None:
            self.reference_union = unary_union(list(self.dict_reference.values()))
        return self.reference_union

    def _postprocess_preresult(self, preresult, geom_thematic):
        """
        Postprocess the preresult with the following actions to create the final result
        *Corrections for areas that differ more than the relevant distance
        *slivers
        *Inner holes (donuts) /multipolygons
        *validity
        *Circles (Polspy-popper)
        *Null/Empty-values

        Args:
            preresult (list): An existing list with all the elements of the preresult
            geom_thematic (BaseGeometry): The input geometry

        Returns:
            tuple: A tuple containing the resulting output geometries:

            *   result (BaseGeometry): The resulting output geometry
            *   result_diff (BaseGeometry): The resulting difference output geometry
            *   result_diff_plus (BaseGeometry): The resulting positive difference output
                geometry
            *   result_diff_min (BaseGeometry): The resulting negative difference output
                geometry
        """
        # Process array
        result = []
        geom_preresult = make_valid(unary_union(preresult))
        geom_thematic = make_valid(geom_thematic)

        # Corrections for areas that differ more than the relevant distance
        geom_thematic_dissolved = buffer_pos(
            buffer_neg(
                buffer_pos(geom_preresult, CORR_DISTANCE),
                2 * CORR_DISTANCE,
            ),
            CORR_DISTANCE,
        )
        # geom_symdiff = self._safe_symmetric_difference(geom_thematic,
        # geom_thematic_dissolved)
        geom_diff_add = safe_difference(geom_thematic, geom_thematic_dissolved)
        geom_diff_delete = safe_difference(geom_thematic_dissolved, geom_thematic)
        geom_diff_removed = safe_difference(
            geom_thematic_dissolved,
            safe_intersection(
                geom_diff_delete,
                buffer_neg_pos(geom_diff_delete, self.buffer_distance()),
            ),
        )
        geom_diff_removed_added = safe_union(
            geom_diff_removed,
            safe_intersection(
                geom_diff_add,
                buffer_neg_pos(geom_diff_add, self.buffer_distance()),
            ),
        )
        geom_thematic_preresult = buffer_pos(
            buffer_neg(
                buffer_pos(geom_diff_removed_added, CORR_DISTANCE),
                2 * CORR_DISTANCE,
            ),
            CORR_DISTANCE,
        )
        # Correction for Inner holes(donuts) / multipolygons
        # fill and remove gaps
        geom_thematic_cleaned_holes = geom_thematic_preresult
        ix_part = 1
        for part in get_parts(geom_thematic_preresult):
            exterior_ring = get_exterior_ring(part)
            exterior_polygon = polygons([exterior_ring])[0]
            empty_buffered_exterior_polygon = buffer_neg(
                exterior_polygon, self.buffer_distance()
            ).is_empty
            if (
                ix_part > 1
                and empty_buffered_exterior_polygon
                and not exterior_polygon.is_empty
            ):
                geom_thematic_cleaned_holes = safe_difference(
                    geom_thematic_cleaned_holes, exterior_polygon
                )
            num_interior_rings = get_num_interior_rings(part)
            if num_interior_rings > 0:
                ix = 0
                while ix < num_interior_rings:
                    interior_ring = get_interior_ring(part, ix)
                    interior_polygon = polygons([interior_ring])[0]

                    empty_buffered_interior_ring = buffer_neg(
                        interior_polygon, self.buffer_distance()
                    ).is_empty
                    if empty_buffered_interior_ring:
                        geom_thematic_cleaned_holes = safe_union(
                            geom_thematic_cleaned_holes, interior_polygon
                        )
                    ix = ix + 1
            ix_part = ix_part + 1

        geom_thematic_result = buffer_pos(
            buffer_neg(
                buffer_pos(geom_thematic_cleaned_holes, CORR_DISTANCE),
                2 * CORR_DISTANCE,
            ),
            CORR_DISTANCE,
        )
        geom_thematic_result = make_valid(remove_repeated_points(geom_thematic_result))
        # Correction for circles

        # calculate ratio to see if it is a circle, and keep the original geometry if a
        # circle: (Polsby-popper score)
        if not (geom_thematic.is_empty or geom_thematic is None):
            if (
                4 * pi * (geom_thematic.area / (geom_thematic.length**2))
                > THRESHOLD_CIRCLE_RATIO
            ):
                self.feedback_warning(
                    "Circle: -->resulting geometry = original geometry"
                )
                geom_thematic_result = geom_thematic
        # Correction for empty preresults
        if geom_thematic_result.is_empty or geom_thematic_result is None:
            self.feedback_warning(
                "Empty result: -->resulting geometry = original geometry"
            )
            geom_thematic_result = geom_thematic

        # group all initial multipolygons into a new resulting dictionary
        result.append(geom_thematic_result)

        # create all resulting geometries
        geom_thematic_result = make_valid(unary_union(result))
        geom_result_diff = buffer_pos(
            buffer_neg(
                safe_symmetric_difference(geom_thematic_result, geom_thematic),
                CORR_DISTANCE,
            ),
            CORR_DISTANCE,
        )
        geom_result_diff_plus = safe_difference(geom_thematic_result, geom_thematic)
        geom_result_diff_min = safe_difference(geom_thematic, geom_thematic_result)
        result_diff_plus = geom_result_diff_plus
        result_diff_min = geom_result_diff_min
        result = geom_thematic_result
        result_diff = geom_result_diff
        return result, result_diff, result_diff_plus, result_diff_min

    @staticmethod
    def add_multi_polygons_from_geom_to_array(geom: BaseGeometry, array):
        """
        Append valid polygons and multipolygons extracted from a given geometry to an
        existing array.

        Args:
            geom (BaseGeometry): The input geometry to process.
            array (list): An existing list to store valid polygons and multipolygons.

        Returns:
            list: A list containing valid polygons and multipolygons extracted from the
                input geometry.
        """
        if geom.is_empty or geom is None:
            # If the input geometry is empty or None, do nothing.
            pass
        else:
            # Create a GeometryCollection from the input geometry.
            geometry_collection = GeometryCollection(geom)
            for g in geometry_collection.geoms:
                # Ensure each sub-geometry is valid.
                g = make_valid(g)
                if str(g.geom_type) in ["Polygon", "MultiPolygon"]:
                    # Append valid polygons and multipolygons to the array.
                    array.append(g)
        return array

    @staticmethod
    def partition(geom, delta):
        """
        Filters a computed grid of partitions (generated by `_grid_bounds`) based on
        intersection with a geometric object (`geom`).

        Args:
            geom (BaseGeometry): The geometric object to check for intersection
                with partitions.
            delta (float): The distance between partitions (same value used in
                `_grid_bounds`).

        Returns:
            list: A filtered list of Polygon objects representing the partitions
                overlapping the original geometric object.
        """
        prepared_geom = prep(geom)
        partitions = grid_bounds(geom, delta)
        filtered_grid = list(filter(prepared_geom.intersects, partitions))
        return filtered_grid
