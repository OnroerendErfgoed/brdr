import hashlib
import json
import logging
import math
import os.path
import uuid
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, date
from decimal import Decimal
from io import BytesIO

import geopandas as gpd
import requests
from geojson import Feature
from geojson import FeatureCollection
from pyproj import CRS
from shapely import GeometryCollection
from shapely import make_valid
from shapely import node
from shapely import polygonize
from shapely.geometry import shape
from shapely.geometry.base import BaseGeometry

from brdr.constants import AREA_ATTRIBUTE
from brdr.constants import DEFAULT_CRS
from brdr.constants import DOWNLOAD_LIMIT
from brdr.constants import PERIMETER_ATTRIBUTE
from brdr.constants import SHAPE_INDEX_ATTRIBUTE
from brdr.constants import STABILITY
from brdr.constants import ZERO_STREAK
from brdr.enums import DiffMetric
from brdr.geometry_utils import buffer_neg, gml_response_to_geojson
from brdr.geometry_utils import buffer_pos
from brdr.geometry_utils import from_crs
from brdr.geometry_utils import get_bbox
from brdr.geometry_utils import get_partitions
from brdr.geometry_utils import get_shape_index
from brdr.geometry_utils import safe_intersection
from brdr.geometry_utils import safe_unary_union
from brdr.geometry_utils import to_crs
from brdr.geometry_utils import total_vertex_distance
from brdr.logger import LOGGER
from brdr.typings import ProcessResult

log = logging.getLogger(__name__)


def _parse_feature_collection_response(
    response, base_url=None, params=None, request_timeout=60
):
    """
    Parse an HTTP response into a GeoJSON FeatureCollection-like dictionary.

    Strategy:
    1. Try JSON first (preferred for OGC API Features and JSON-capable WFS).
    2. Fallback to in-memory GML/XML parsing using GeoPandas.
    3. Last fallback: re-fetch through `gml_response_to_geojson` when URL/params are available.
    """
    try:
        data = response.json()
        if isinstance(data, dict):
            return data
    except Exception:
        pass

    try:
        gdf = gpd.read_file(BytesIO(response.content))
        return json.loads(gdf.to_json())
    except Exception:
        if base_url is not None:
            return gml_response_to_geojson(
                base_url, params or {}, timeout=request_timeout
            )
        raise


def _request_with_outputformat_fallback(
    base_url, params=None, headers=None, request_timeout=60
):
    """
    Request helper that retries once without `outputFormat` when unsupported.
    """
    local_params = (params or {}).copy()
    response = requests.get(
        base_url, params=local_params, headers=headers, timeout=request_timeout
    )
    if response.ok:
        return response, local_params

    if "outputFormat" in local_params:
        LOGGER.info(
            "Request failed with outputFormat=%s; retrying without outputFormat.",
            local_params.get("outputFormat"),
        )
        local_params.pop("outputFormat", None)
        response_retry = requests.get(
            base_url, params=local_params, headers=headers, timeout=request_timeout
        )
        if response_retry.ok:
            return response_retry, local_params
        response_retry.raise_for_status()

    response.raise_for_status()
    return response, local_params


def get_geojsons_from_process_results(
    process_results: dict[str | int, dict[float, ProcessResult]],
    crs: CRS,
    id_field: str,
    series_prop_dict: dict[str | int, dict[float, str | int]] = None,
    geom_attributes=True,
):
    """
    Convert a series of process results to a GeoJSON feature collection.

    Args:
        process_results (dict): Dictionary containing process results.
        crs (str): Coordinate reference system.
        id_field (str): Field name for the ID.
        series_prop_dict (dict, optional): Dictionary containing series properties.
        geom_attributes (bool, optional): Whether to include geometry attributes.

    Returns:
        dict: Dictionary of GeoJSON feature collections.
    """
    features_list_dict = {}

    for theme_id, results_dict in process_results.items():

        prop_dict = dict(series_prop_dict or {}).get(theme_id, {})
        for relevant_distance, process_result in results_dict.items():
            properties = prop_dict.get(relevant_distance, {})
            properties[id_field] = theme_id
            properties.update(process_result["properties"])
            for results_type, value in process_result.items():
                if not isinstance(value, BaseGeometry):
                    continue
                if results_type not in features_list_dict:
                    features_list_dict[results_type] = []

                feature = _feature_from_geom(
                    value, theme_id, properties, geom_attributes
                )
                features_list_dict[results_type].append(feature)

    crs_geojson = {"type": "name", "properties": {"name": from_crs(crs)}}
    geojsons = {
        result_type: FeatureCollection(features, crs=crs_geojson)
        for result_type, features in features_list_dict.items()
    }
    return geojsons


def deep_merge(dict_a, dict_b):
    for key, value in dict_b.items():
        if key in dict_a and isinstance(dict_a[key], dict) and isinstance(value, dict):
            deep_merge(dict_a[key], value)
        else:
            dict_a[key] = value
    return dict_a


def _feature_from_geom(
    geom: BaseGeometry,
    feature_id: str | int,
    properties: dict = None,
    geom_attributes=True,
) -> Feature:
    """
    Convert a geometry to a GeoJSON feature.

    Args:
        geom (BaseGeometry): The geometry to convert.
        feature_id: a (unique) id for the feature
        properties (dict, optional): The properties to include in the feature.
        geom_attributes (bool, optional): Whether to include geometry attributes.

    Returns:
        Feature: The GeoJSON feature.
    """
    properties = dict(properties or {})
    if geom_attributes:
        area = geom.area
        perimeter = geom.length
        properties[AREA_ATTRIBUTE] = area
        properties[PERIMETER_ATTRIBUTE] = perimeter
        properties[SHAPE_INDEX_ATTRIBUTE] = get_shape_index(area, perimeter)
    return Feature(geometry=geom, id=feature_id, properties=properties)


def write_featurecollection_to_geopackage(
    output_path, featurecollection, layer_name="resulting_layer"
):
    # 1. Change to GeoDataFrame
    crs_name = featurecollection.get("crs", {}).get("properties", {}).get("name")
    gdf = gpd.GeoDataFrame.from_features(featurecollection, crs=crs_name)

    # 2. Set CRS
    if gdf.crs is None:
        gdf.set_crs(DEFAULT_CRS, inplace=True)

    # 3. serialize dicts and lists to string
    for col in gdf.columns:
        if col == "geometry":
            continue

        sample = gdf[col].dropna().iloc[0] if not gdf[col].dropna().empty else None
        if isinstance(sample, (dict, list)):
            gdf[col] = gdf[col].apply(
                lambda x: json.dumps(x, ensure_ascii=False) if x is not None else None
            )

    # Write to geopackage
    folder = os.path.dirname(output_path)
    if folder and not os.path.exists(folder):
        os.makedirs(folder, exist_ok=True)
    gdf.to_file(
        output_path,
        driver="GPKG",
        layer=layer_name,
        # mode="a",
        engine="pyogrio",
        # ,overwrite_layer=True
    )


def geojson_serializor(obj):
    if isinstance(obj, (datetime, date)):
        return obj.isoformat()
    if isinstance(obj, Decimal):
        return float(obj)
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")


def write_featurecollection_to_geojson(output_path, featurecollection):
    features = featurecollection.get("features", [])

    # Determine all columns who need a json.dumps (lists and dicts)
    keys_to_stringify = set()
    if features:
        sample_props = features[0]["properties"]
        for k, v in sample_props.items():
            if isinstance(v, (dict, list)):
                keys_to_stringify.add(k)

    # Serializes these columns
    for feature in features:
        props = feature["properties"]
        for key in keys_to_stringify:
            if props[key] is not None:
                props[key] = json.dumps(props[key], ensure_ascii=False)

    # create/write file
    parent = os.path.dirname(output_path)
    os.makedirs(parent, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(
            featurecollection,
            f,
            indent=2,
            default=geojson_serializor,
            ensure_ascii=False,
        )


def flatten_iter(lst):
    while any(isinstance(i, list) for i in lst):
        lst = [
            item
            for sublist in lst
            for item in (sublist if isinstance(sublist, list) else [sublist])
        ]
    return lst


def polygonize_reference_data(dict_ref):
    """
    Create a new dictionary with non-overlapping polygons based on a reference data dictionary.

    Args:
        dict_ref (dict): Dictionary of reference geometries.

    Returns:
        dict: Dictionary of non-overlapping polygons.
    """
    arr_ref = []
    for key in dict_ref:
        arr_ref.append((make_valid(dict_ref[key])))
    arr_ref = node(GeometryCollection(arr_ref))
    polygonized = polygonize(arr_ref.geoms)
    dict_ref = {}
    i = 0
    for g in polygonized.geoms:
        if str(g.geom_type) in ["Polygon", "MultiPolygon"]:
            i = i + 1
            key = str(i)  # unique keys will be added (reference_id is lost)
            dict_ref[key] = make_valid(g)
        else:
            logging.warning("geom excluded: " + str(g))
    return dict_ref


def coverage_ratio(values, min_val=0, bin_count=10):
    max_val = max(values)
    if len(values) == 0 or max_val <= 0:
        return 0.0
    bin_size = round((max_val - min_val) / bin_count, 2)
    bins = np.arange(min_val, max_val + bin_size, bin_size)

    # Tel hoeveel bins minstens één waarde bevatten
    bin_counts, _ = np.histogram(values, bins)
    filled_bins = np.count_nonzero(bin_counts)

    total_bins = len(bins) - 1  # histogram geeft n-1 bins
    return filled_bins / total_bins


import numpy as np


def determine_stability2(
    x_coords, y_coords, max_derivative_threshold=100, num_plateaus=10
):
    """
    Determine the stability and indicators ('zero_streaks') based on the derivative of the difference measurement (y) of a xy-plot

    Args:
        x_coords (numpy.ndarray): The x values of the plot: relevant distances
        y_coords (numpy.ndarray): The y values of the plot: diffence measurement

    Returns:
        tuple: A tuple containing:
            - list: List of breakpoints (extremes).
            - list: List of zero_streaks.
    """

    # Calculate derivatives (slopes)
    dx = np.diff(x_coords)
    dy = np.diff(y_coords)
    slopes = dy / dx

    found_segments = []

    # 1. Search for the best plateaus
    for i in range(1, len(slopes)):
        slope_before = slopes[i - 1]
        best_score_for_i = -1
        best_seg_for_i = None

        # Scan forward to find the optimal end point for a plateau starting at x[i]
        for j in range(i, len(slopes)):
            avg_slope_seg = (y_coords[j + 1] - y_coords[i]) / (
                x_coords[j + 1] - x_coords[i]
            )
            length_seg = x_coords[j + 1] - x_coords[i]

            # Condition: Average slope must stay below threshold and show significant deceleration
            if (
                avg_slope_seg > max_derivative_threshold
                or avg_slope_seg > slope_before * 0.7
            ):
                break

            # Score Calculation (Max 100)
            # Flatness (50pts): Derivative 0 = 50pts, Derivative at threshold = 0pts
            flatness_score = max(
                0,
                (max_derivative_threshold - avg_slope_seg)
                / max_derivative_threshold
                * 50,
            )

            # Length (50pts): Normalized to 1.0 unit of x for full score
            length_score = min(50, (length_seg / 1.0) * 50)

            total_score = flatness_score + length_score

            if total_score > best_score_for_i:
                best_score_for_i = total_score
                best_seg_for_i = (x_coords[i], x_coords[j + 1], total_score)

        if best_seg_for_i:
            found_segments.append(best_seg_for_i)

    # Filter overlapping segments and keep only the top ones
    found_segments = sorted(found_segments, key=lambda k: k[2], reverse=True)
    top_plateaus = {}
    chosen_areas = []

    for start, stop, score in found_segments:
        # Avoid overlapping: a new plateau cannot overlap with a better one already chosen
        if not any(not (stop <= u[0] or start >= u[1]) for u in chosen_areas):
            top_plateaus[start] = (start, stop, score)
            chosen_areas.append((start, stop))
            if len(top_plateaus) >= num_plateaus:
                break

    # 2. Build the final dictionary for ALL x-values
    output = {}
    for val_x in x_coords:
        is_stability_point = val_x in top_plateaus

        if is_stability_point:
            start, stop, score = top_plateaus[val_x]
            center_value = round(start + (stop - start) / 2, 2)
            streak_tuple = (
                round(start, 2),  # start_value
                round(stop, 2),  # end_value
                center_value,  # center_center value
                round(score, 2),  # score
            )
        else:
            streak_tuple = None

        output[val_x] = {
            "brdr_stability": is_stability_point,
            "brdr_zero_streak": streak_tuple,
        }

    return output


def determine_stability(x, y):
    """
    Determine the stability and indicators ('zero_streaks') based on the derivative of the difference measurement (y) of a xy-plot

    Args:
        x (numpy.ndarray): The x values of the plot: relevant distances
        y (numpy.ndarray): The y values of the plot: diffence measurement

    Returns:
        tuple: A tuple containing:
            - list: List of breakpoints (extremes).
            - list: List of zero_streaks.
    """
    derivative = _numerical_derivative(x, y)
    max_zero_streak_score = x[-1] - x[0]
    start_streak = None
    streak = 0
    write_zero_streak = False
    dict_stability = dict()
    for i in range(1, len(x)):
        dict_stability[x[i - 1]] = dict()
        dict_stability[x[i - 1]][STABILITY] = False
        dict_stability[x[i - 1]][ZERO_STREAK] = None
        if round(derivative[i], 2) == 0:  # noqa
            dict_stability[x[i - 1]][STABILITY] = True
            streak = streak + 1
            if start_streak is None:
                start_streak = x[i - 1]
        elif streak != 0:
            write_zero_streak = True
        if start_streak is not None and (write_zero_streak or len(x) - 1 == i):
            end_streak = x[i - 1]
            if derivative[i] == 0:
                end_streak = x[i]
            score_streak = round(
                (end_streak - start_streak) * 100 / max_zero_streak_score, 2
            )
            center_streak = start_streak + (end_streak - start_streak) / 2

            dict_stability[start_streak][ZERO_STREAK] = (
                start_streak,
                end_streak,
                center_streak,
                score_streak,
            )

            streak = 0
            start_streak = None
            write_zero_streak = False
            logging.debug("end_streak")

    return dict_stability


def _numerical_derivative(x, y):
    """
    Calculate the numerical derivative of a graph.

    Parameters:
      x (numpy.ndarray): The x values of the graph.
      y (numpy.ndarray): The y values of the graph.

    Returns:
      numpy.ndarray: The derivative of y with respect to x.
    """

    dx = np.diff(x)
    dy = np.diff(y)
    derivative = dy / dx
    derivative = np.insert(derivative, 0, 0)
    # derivative = np.append(derivative, 0)

    return derivative


def get_geometry_difference_metrics_from_processresults(
    dict_processresult: dict[float, ProcessResult],
    geom_thematic: BaseGeometry,
    reference_union: BaseGeometry,
    diff_metric: DiffMetric = DiffMetric.SYMMETRICAL_AREA_CHANGE,
):
    """
    Calculates a dictionary containing difference metrics for thematic geometry based on a distance series.

    Parameters:
    dict_processresult (dict): A dictionary where keys are thematic IDs and values are dictionaries mapping relevant distances to ProcessResult objects.
    geom_thematic (BaseGeometry): A dictionary where keys are thematic IDs and values are BaseGeometry objects representing the original geometries.
    reference_union (BaseGeometry): unioned reference geometries
    diff_metric (DiffMetric, optional): The metric to use for calculating differences. Default is DiffMetric.CHANGES_AREA.

    Returns:
    dict: A dictionary where keys are relevant distances to calculated difference metrics.
    """

    diffs = {}

    for rel_dist in dict_processresult:
        processresult = dict_processresult.get(rel_dist, {})
        if geom_thematic.geom_type in (
            "LineString",
            "MultiLineString",
        ):
            diff_metric = DiffMetric.LENGTH_ADDED_AND_REMOVED
        elif geom_thematic.geom_type in (
            "Point",
            "MultiPoint",
        ):
            diff_metric = DiffMetric.TOTAL_DISTANCE
        diff = get_geometry_difference_metrics_from_processresult(
            processresult, geom_thematic, reference_union, diff_metric
        )

        diffs[rel_dist] = diff
    return diffs


def get_geometry_difference_metrics_from_processresult(
    processresult, geom_thematic, reference_union, diff_metric
):
    result = processresult.get("result")
    result_diff = processresult.get("result_diff")
    result_diff_min = processresult.get("result_diff_min")
    result_diff_plus = processresult.get("result_diff_plus")
    diff = 0
    original = geom_thematic
    if result_diff is None or result_diff.is_empty or result is None or result.is_empty:
        diff = 0
    elif diff_metric == DiffMetric.SYMMETRICAL_AREA_CHANGE:
        diff = result_diff.area
    elif diff_metric == DiffMetric.SYMMETRICAL_AREA_PERCENTAGE_CHANGE:
        diff = result_diff.area
        if original.area != 0:
            diff = diff * 100 / original.area
    elif diff_metric == DiffMetric.AREA_CHANGE:
        diff = result.area - original.area
    elif diff_metric == DiffMetric.AREA_PERCENTAGE_CHANGE:
        diff = result.area - original.area
        if original.area != 0:
            diff = diff * 100 / original.area
    elif diff_metric == DiffMetric.LENGTH_CHANGE:
        diff = result.length - original.length
    elif diff_metric == DiffMetric.LENGTH_PERCENTAGE_CHANGE:
        diff = result.length - original.length
        if original.length != 0:
            diff = diff * 100 / original.length
    elif diff_metric == DiffMetric.LENGTH_REMOVED:
        if not result_diff_min is None or result_diff_min.is_empty:
            diff = result_diff_min.length
    elif diff_metric == DiffMetric.LENGTH_ADDED_AND_REMOVED:
        if not result_diff_min is None or result_diff_min.is_empty:
            diff_min = result_diff_min.length
        else:
            diff_min = 0
        if not result_diff_plus is None or result_diff_plus.is_empty:
            diff_plus = result_diff_plus.length
        else:
            diff_plus = 0
        diff = diff_plus + diff_min
    elif diff_metric == DiffMetric.REFERENCE_USAGE:
        if not reference_union is None and not reference_union.is_empty:
            reference_union_buffer = buffer_pos(reference_union, 0.01)
            result_buffer = buffer_pos(result, 0.01)
            reference_usage_geom = safe_intersection(
                result_buffer, reference_union_buffer
            )
        else:
            reference_usage_geom = None
        if reference_usage_geom is not None and not reference_usage_geom.is_empty:
            diff = reference_usage_geom.area
        else:
            diff = 0
    elif diff_metric == DiffMetric.TOTAL_DISTANCE:
        diff = total_vertex_distance(original, result, bidirectional=False)

    # round, so the detected changes are within 10cm, 10cm² or 0.1%
    diff = round(abs(diff), 1)
    return diff


def fetch_all_ogc_features(
    base_url,
    params=None,
    headers=None,
    max_pages=math.inf,
    max_workers=None,
    request_timeout=60,
):
    """
    Fetches all features from an OGC Feature API using parallel requests where possible.

    If the API provides 'numberMatched', the function calculates all offsets and
    fetches pages in parallel. If not, it falls back to sequential cursor-based
    pagination via 'next' links.

    Parameters
    ----------
    base_url : str
        URL of the /items endpoint of the OGC API.
    params : dict, optional
        Query parameters for the API request. Default is None.
    headers : dict, optional
        Optional headers (e.g., Accept: application/json).
    max_pages : int or math.inf, optional
        Maximum number of pages to fetch. Default is math.inf.
    max_workers : int, optional
        Maximum number of concurrent threads for parallel fetching. Default is None.

    Returns
    -------
    list
        A list containing all features retrieved from the API.
    """
    if params is None:
        params = {}
    if headers is None:
        headers = {"Accept": "application/json"}

    requested_limit = params.get("limit", DOWNLOAD_LIMIT)

    # 1. Initial request to determine total count and pagination type
    response, params = _request_with_outputformat_fallback(
        base_url, params=params, headers=headers, request_timeout=request_timeout
    )
    data = _parse_feature_collection_response(
        response,
        base_url=base_url,
        params=params,
        request_timeout=request_timeout,
    )

    all_features = list(data.get("features", []))
    total_matched = data.get("numberMatched")
    number_returned = data.get("numberReturned", len(data.get("features", [])))
    # Use the actual returned page size as effective step for pagination when the
    # server silently enforces a lower max feature count than the requested limit.
    if number_returned and number_returned > 0:
        effective_limit = number_returned
    else:
        effective_limit = requested_limit

    # Check for 'next' link for cursor-based fallback
    next_link = next(
        (link["href"] for link in data.get("links", []) if link["rel"] == "next"),
        None,
    )

    # 2. Parallel Path: If we know the total and it's more than one page
    if total_matched and total_matched > number_returned and max_pages > 0:
        LOGGER.info(
            f"Total features: {total_matched}. Switching to parallel offset-based retrieval."
        )

        # Calculate offsets
        offsets = list(range(number_returned, total_matched, effective_limit))
        # Limit offsets by max_pages only when finite
        if isinstance(max_pages, (int, float)) and math.isfinite(max_pages):
            offsets = offsets[: int(max_pages)]

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_offset = {}
            for offset in offsets:
                local_params = params.copy()
                local_params["startIndex"] = offset
                local_params["limit"] = effective_limit
                future_to_offset[
                    executor.submit(
                        _request_with_outputformat_fallback,
                        base_url,
                        local_params,
                        headers,
                        request_timeout,
                    )
                ] = offset

            for future in as_completed(future_to_offset):
                try:
                    resp, used_params = future.result()
                    page_data = _parse_feature_collection_response(
                        resp,
                        base_url=base_url,
                        params=used_params,
                        request_timeout=request_timeout,
                    )
                    all_features.extend(page_data.get("features", []))
                except Exception as exc:
                    LOGGER.error(
                        f"Offset {future_to_offset[future]} generated an exception: {exc}"
                    )

        return all_features

    # 3. Sequential Path: Fallback for Cursor-based pagination
    elif next_link and max_pages > 0:
        LOGGER.info(
            "Parallel retrieval not possible (no total count). Falling back to sequential cursor pagination."
        )
        url = next_link
        page = 1
        seen_next_urls = set()
        while url and page <= max_pages:
            if url in seen_next_urls:
                LOGGER.warning(
                    "Detected cyclic 'next' pagination URL. Stopping sequential pagination."
                )
                break
            seen_next_urls.add(url)
            LOGGER.debug(f"Fetching page {page}: {url}")
            # For next_links, parameters are usually already in the URL
            response = requests.get(url, headers=headers, timeout=request_timeout)
            response.raise_for_status()
            data = _parse_feature_collection_response(
                response,
                base_url=url,
                params={},
                request_timeout=request_timeout,
            )

            all_features.extend(data.get("features", []))

            url = next(
                (
                    link["href"]
                    for link in data.get("links", [])
                    if link["rel"] == "next"
                ),
                None,
            )
            page += 1

    return all_features


def make_feature_collection(features):
    """
    Create a GeoJSON FeatureCollection from a list of GeoJSON features.

    Parameters
    ----------
    features : list[dict]
        List of GeoJSON features (`{"type": "Feature", ...}`).

    Returns
    -------
    dict
        GeoJSON FeatureCollection dictionary.
    """
    return {"type": "FeatureCollection", "features": features}


def get_collection(
    url, params=None, max_pages=math.inf, max_workers=None, request_timeout=60
):
    """
    Fetches a collection of features from a (paginated) API endpoint (OGC Feature API or WFS).

    This function retrieves a collection of features from a URL that supports
    pagination using a `startIndex` parameter. It iteratively retrieves features in
    chunks of the specified `limit` until no more features are available.

    Args:
        url (str): The base URL of the API endpoint.
        params: parameters for the request

    Returns:
        dict: A dictionary representing the complete GeoJSON feature collection.
            This might be truncated if the total number of features exceeds the
            limitations of the API or server.
    """
    if params is None:
        params = {}
    features = fetch_all_ogc_features(
        base_url=url,
        params=params,
        max_pages=max_pages,
        max_workers=max_workers,
        request_timeout=request_timeout,
    )
    return make_feature_collection(features)


def geojson_geometry_to_shapely(geojson_geometry):
    """
    Convert a GeoJSON geometry object to a Shapely geometry.

    Parameters
    ----------
    geojson_geometry : dict
        GeoJSON geometry object.

    Returns
    -------
    BaseGeometry
        Converted Shapely geometry.
    """
    return shape(geojson_geometry)


def geojson_to_dicts(collection, id_property=None):
    """
    Converts a GeoJSON collection into dictionaries of geometries and properties.

    Parameters:
    collection (dict): The GeoJSON collection to convert.
    id_property (str): The property name to use as the key for the dictionaries. if None, the feature-identifier 'id' is used when available.

    Returns:
    tuple: Two dictionaries:
        - data_dict (dict): A dictionary where keys are the id_property values and values are the geometries.
        - data_dict_properties (dict): A dictionary where keys are the id_property values and values are the properties.
    """
    data_dict = {}
    data_dict_properties = {}
    if collection is None or "features" not in collection:
        return data_dict, data_dict_properties
    for f in collection["features"]:
        if not id_property:
            if "id" in f:
                key = f["id"]
            else:
                raise KeyError(
                    "Feature-identifier 'id' not found in GeoJson FeatureCollection. Please provide a Geojson with feature-identifiers or define id_property (name of attribute with unique ids) from the properties"
                )
        else:
            if "properties" in f and id_property in f["properties"]:
                key = f["properties"][id_property]
            else:
                raise KeyError(
                    f"Identifier '{id_property}' not found in properties of GeoJson FeatureCollection"
                )
        geom = shape(f["geometry"])
        data_dict[key] = make_valid(geom)
        data_dict_properties[key] = f["properties"]
    return data_dict, data_dict_properties


def get_collection_by_partition(
    url,
    params,
    geometry,
    partition=1000,
    crs=DEFAULT_CRS,
    max_workers=None,
    max_pages=math.inf,
    request_timeout=60,
):
    """
    Retrieves a collection of geographic data by partitioning the input geometry.

    This function parallelizes data retrieval by splitting the geometry into
    smaller parts and fetching each partition concurrently using threads.

    Parameters
    ----------
    url : str
        The base URL for the OGC API Features data source.
    params : dict
        The query parameters for the API request.
    geometry : object
        The geometric area (Shapely object) to partition and retrieve data for.
        If None or empty, retrieves data for the entire area.
    partition : int, optional
        The number of partitions to divide the geometry into. Default is 1000.
        If less than 1, no partitioning is performed.
    crs : str, optional
        The coordinate reference system to use. Default is DEFAULT_CRS.
    max_workers : int, optional
        The maximum number of concurrent threads to use for parallel requests.
        Default is None.

    Returns
    -------
    dict
        A GeoJSON-like collection of geographic data containing merged features
        from all partitions.
    """

    crs = to_crs(crs)
    collection = {}

    if geometry is None or geometry.is_empty:
        return get_collection(
            url=url,
            params=params,
            max_pages=max_pages,
            max_workers=max_workers,
            request_timeout=request_timeout,
        )

    if partition < 1:
        local_params = params.copy()
        local_params["bbox"] = get_bbox(geometry)
        local_params["bbox-crs"] = from_crs(crs)
        return get_collection(
            url=url,
            params=local_params,
            max_pages=max_pages,
            max_workers=max_workers,
            request_timeout=request_timeout,
        )

    # Prepare partitions
    geoms = get_partitions(geometry, partition)
    features_list = []

    # Use ThreadPoolExecutor for parallel I/O-bound API requests
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_geom = {}

        for g in geoms:
            # Create a shallow copy of params to avoid race conditions
            thread_params = params.copy()
            thread_params["bbox"] = get_bbox(g)
            thread_params["bbox-crs"] = from_crs(crs)

            # Submit the request to the thread pool
            future = executor.submit(
                get_collection,
                url=url,
                params=thread_params,
                max_pages=max_pages,
                max_workers=max_workers,
                request_timeout=request_timeout,
            )
            future_to_geom[future] = g

        # Collect results as they complete
        for future in as_completed(future_to_geom):
            try:
                coll = future.result()
                if coll and "features" in coll:
                    # If this is the first valid response, use it as the template for the collection
                    if not collection:
                        collection = dict(coll)
                        collection["features"] = []

                    features_list.extend(coll["features"])
            except Exception as exc:
                # Retrieve the specific geometry that caused the failure for debugging
                geom_info = future_to_geom[future]
                print(f"Partition request generated an exception: {exc}")
                raise exc

    # Re-attach all collected features to the template
    if collection:
        unique_features = deduplicate_features(features_list)
        collection["features"] = unique_features
    # TODO; only return features that intersect with the 'geometry'?

    return collection


def deduplicate_features(features):
    """
    Removes duplicate features using ID or a content-based hash as fallback.

    Parameters
    ----------
    features : list
        List of GeoJSON-like feature dictionaries.

    Returns
    -------
    list
        Deduplicated list of features.
    """
    seen_hashes = set()
    unique_features = []

    for feat in features:
        # 1. Try to get the official ID
        feat_id = feat.get("id")

        # 2. If no ID, create a hash of the content
        if feat_id is None:
            # We use a stable JSON string to represent the feature content
            content = json.dumps(
                {"g": feat.get("geometry"), "p": feat.get("properties")}, sort_keys=True
            )
            feat_id = hashlib.md5(content.encode()).hexdigest()

        if feat_id not in seen_hashes:
            unique_features.append(feat)
            seen_hashes.add(feat_id)

    return unique_features


def build_reverse_index_wkb(d: dict):
    return {g.wkb: k for k, g in d.items()}


def is_brdr_observation(brdr_observation):
    """
    Check whether a value matches the expected BRDR observation structure.

    Parameters
    ----------
    brdr_observation : Any
        Candidate observation object.

    Returns
    -------
    bool
        `True` if structure is valid, else `False`.
    """
    if brdr_observation is None or not isinstance(brdr_observation, dict):
        return False
    if brdr_observation.keys() >= {
        "alignment_date",
        "brdr_version",
        "reference_source",
        "full",
        "area",
        "reference_features",
        "reference_od",
    }:
        return True
    return False


def union_process_result(process_result: ProcessResult) -> ProcessResult:
    """
    Perform a unary union on all geometries within a ProcessResult object.

    This function iterates through all annotated geometric fields in the
    ProcessResult, attempts to merge overlapping or adjacent geometries
    into a single unified structure using a safe unary union operation,
    and updates the result in place.

    Parameters
    ----------
    process_result : ProcessResult
        A dictionary-like or dataclass object containing geometric data
        (typically resulting from a polygon processing operation).
        The keys correspond to specific geometric categories or layers.

    Returns
    -------
    ProcessResult
        The updated ProcessResult object where each valid, non-empty
        geometry has been replaced by its unioned representation.

    Notes
    -----
    - Empty geometries or those that are not instances of `BaseGeometry`
      are ignored.
    - The function utilizes `safe_unary_union` to handle potential
      topological inconsistencies or invalid geometries that might
      cause a standard union to fail.
    - Fields are accessed based on the type annotations of the
      `ProcessResult` class.

    See Also
    --------
    safe_unary_union : The underlying function used for geometric merging.
    ProcessorConfig : Settings that define how these results are generated.
    """
    for key in ProcessResult.__annotations__:
        if key in ["properties", "metadata", "observations"]:
            process_result[key] = process_result.get(key, {})  # noqa
            continue
        value = process_result.get(key, GeometryCollection())  # noqa
        if isinstance(value, BaseGeometry):
            value = safe_unary_union(value)
        process_result[key] = value  # noqa

    return process_result


def get_relevant_polygons_from_geom(
    geometry: BaseGeometry, buffer_distance: float, mitre_limit
):
    """
    Get only the relevant parts (polygon) from a geometry.
    Points, Lines and Polygons smaller than relevant distance are excluded from the
    result
    """
    if not geometry or geometry.is_empty:
        # If the input geometry is empty or None, do nothing.
        return geometry
    else:
        geometry = safe_unary_union(geometry)
        # Create a GeometryCollection from the input geometry.
        geometry_collection = GeometryCollection(geometry)
        array = []
        for g in geometry_collection.geoms:
            # Ensure each sub-geometry is valid.
            g = make_valid(g)
            if str(g.geom_type) in ["Polygon", "MultiPolygon"]:
                relevant_geom = buffer_neg(
                    g,
                    buffer_distance,
                    mitre_limit=mitre_limit,
                )
                if relevant_geom is not None and not relevant_geom.is_empty:
                    array.append(g)
    return safe_unary_union(array)


def urn_from_geom(geom: BaseGeometry):
    return uuid.UUID(hex=hashlib.sha256(geom.wkb).hexdigest()[::2]).urn
