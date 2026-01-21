from datetime import datetime
from math import inf
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely import wkt, GeometryCollection

from brdr.aligner import Aligner
from brdr.be.grb.enums import GRBType
from brdr.be.grb.loader import GRBActualLoader
from brdr.configs import AlignerConfig
from brdr.constants import RELEVANT_DISTANCE_DECIMALS, DEFAULT_CRS
from brdr.enums import AlignerResultType
from brdr.geometry_utils import (
    buffer_neg,
    geom_to_wkt,
    safe_unary_union,
    buffer_pos,
    to_crs,
    from_crs,
)
from brdr.loader import DictLoader
from brdr.utils import geojson_geometry_to_shapely
from brdr.viz import export_boxplot, export_histogram

METRICS_TXT = "metrics.txt"
STATS_TXT = "stats.txt"

FALSE_POSITIVE_WKT = "falsepositive_wkt"
DOUBT_WKT = "doubt_wkt"
BRDR_WKT = "brdr_wkt"
ORIGINAL_WKT = "original_wkt"

FP_ANALYSIS_CSV_NAME = "fp_analysis.csv"
FP_HISTOGRAM_NAME = "fp_histogram.png"
FP_BOXPLOT_NAME = "fp_boxplot.png"
FP_ESIMATION_COLUMN_NAME = "fp_estimation"


def get_parcel_lists(aligner, geom):
    parcel_100perc_list = []
    observation = aligner.compare_to_reference(geom)
    reference_features_list = []
    if "reference_features" in observation:
        reference_features_list = observation["reference_features"].keys()
        for k, v in observation["reference_features"].items():
            if v["full"]:
                parcel_100perc_list.append(k)
    return list(reference_features_list), parcel_100perc_list


def get_buffer_data(aligner, geom, buffers=[0.1, 0.2, 0.5, 1, 2, 3, 4, 5]):
    buffer_dict = {}
    for buffer in buffers:
        buffer_dict[buffer] = {}
        geom_buffered = buffer_neg(geom, buffer)
        parcel_list, parcel_100perc_list = get_parcel_lists(aligner, geom_buffered)
        buffer_dict[buffer]["geometry"] = geom_buffered
        buffer_dict[buffer]["parcels"] = parcel_list
        buffer_dict[buffer]["parcels_full"] = parcel_100perc_list

    return buffer_dict


def get_coverage_data(aligner, geom, percentages=[0, 1, 5, 10, 50, 90, 95, 99, 100]):

    coverage_dict = {item: {"parcels": [], "count": 0} for item in percentages}

    observation = aligner.compare_to_reference(geom)
    if "reference_features" in observation:
        for k, v in observation["reference_features"].items():
            area_percentage = v["percentage"]
            for p in percentages:
                if area_percentage >= p:
                    coverage_dict[p]["parcels"].append(k)
                    coverage_dict[p]["count"] += 1

    return coverage_dict


def get_false_positive_grb_parcels_dataframe(
    data,
    id_name,
    area_limit=inf,
    processor=None,
    conn_str=None,
    tablename="adp",
    ref_id_name="capakey",
    geometry_name="geom",
    crs=DEFAULT_CRS,
):
    """
    Processes features to find false positives and logs processing metrics.
    """

    # 1. Initialize counters
    starttime = datetime.now()
    metrics = {
        "area_limit": area_limit,
        "total_count": 0,
        "success_ids": [],
        "skipped_ids": [],
        "failed_ids": [],
        "starttime": starttime,
        "endtime": starttime,
    }
    # if not area_limit==inf:
    #     metrics ["area_limit"] = area_limit

    rds = [
        round(k, RELEVANT_DISTANCE_DECIMALS)
        for k in np.arange(0, 510, 10, dtype=int) / 100
    ]
    srid = from_crs(to_crs(crs), format="id")
    dict_coverage_list = {}
    dict_coverage_range = {}
    dict_coverage_0 = {}
    dict_coverage_50 = {}
    dict_buffer_list = {}
    dict_buffer_vip = {}
    dict_buffer_2 = {}
    dict_buffer_5 = {}
    dict_prediction_rd = {}
    dict_prediction = {}
    dict_prediction_wkt = {}
    dict_prediction_state = {}
    dict_prediction_score = {}
    dict_fp = {}
    dict_fp_parcels = {}
    dict_fp_wkt = {}
    dict_doubt_parcels = {}
    dict_doubt_wkt = {}
    dict_original_wkt = {}

    features = data.get("features", []) if data else []
    metrics["total_count"] = len(features)
    print(f"Nr of features: {metrics['total_count']}")
    count = 0

    for f in features:
        count += 1
        print(f"{datetime.now()}: Processing feature {count}/{metrics['total_count']}")
        try:
            key = f["properties"][id_name]
            geom = geojson_geometry_to_shapely(f["geometry"])
            if geom.area > area_limit:
                print(f"Feature {key} skipped based on area {geom.area} m²")
                metrics["skipped_ids"].append(key)
                continue
            print(f"Feature - key: {str(key)}")
            dict_theme = {key: geom}

            aligner_config = AlignerConfig
            aligner = Aligner(crs=crs, processor=processor, config=aligner_config)
            loader = DictLoader(dict_theme)
            aligner.load_thematic_data(loader)
            if not conn_str:
                loader = GRBActualLoader(
                    grb_type=GRBType.ADP, partition=1000, aligner=aligner
                )
            else:
                # SQL-query to get geometry from postgis db
                wkt = buffer_pos(aligner.thematic_data.union, 10)
                sql = f"SELECT {ref_id_name},{geometry_name} from {tablename} where ST_Intersects({geometry_name},ST_SetSRID(ST_GeomFromText('{wkt}'), {str(srid)}))"
                # Read SQL-result into GeoDataFrame
                gdf = gpd.read_postgis(sql, conn_str, geom_col=geometry_name)
                dict_reference = dict(zip(gdf[ref_id_name], gdf[geometry_name]))
                loader = DictLoader(data_dict=dict_reference)
                aligner.load_reference_data(loader)
            aligner.load_reference_data(loader)
            dict_original_wkt = geom.wkt
            # COVERAGE ANALYSIS
            coverage_data = get_coverage_data(
                aligner, geom, percentages=[0, 1, 5, 10, 50, 90, 95, 99, 100]
            )
            parcel_coverage_counts = [
                coverage_data[b]["count"] for b in coverage_data.keys()
            ]
            dict_coverage_list[key] = parcel_coverage_counts
            dict_coverage_range[key] = (
                parcel_coverage_counts[0] - parcel_coverage_counts[-1]
            )
            dict_coverage_0[key] = parcel_coverage_counts[0]
            dict_coverage_50[key] = parcel_coverage_counts[4]

            # BUFFER ANALYSIS
            buffer_data = get_buffer_data(
                aligner, geom, buffers=[0.1, 0.2, 0.5, 1, 2, 3, 4, 5]
            )
            parcel_buffer_counts = [
                len(buffer_data[b]["parcels"]) for b in buffer_data.keys()
            ]
            dict_buffer_list[key] = parcel_buffer_counts
            dict_buffer_vip[key] = parcel_buffer_counts[1]
            # parcels_vip = buffer_data.get(0.2)["parcels"]
            dict_buffer_2[key] = parcel_buffer_counts[4]
            dict_buffer_5[key] = parcel_buffer_counts[-1]

            # BRDR ANALYSIS
            # This is done for every feature seperately
            aligner_result = aligner.evaluate(
                relevant_distances=rds, max_predictions=1, multi_to_best_prediction=True
            )
            dict_predictions_evaluated = aligner_result.get_results(
                aligner=aligner, result_type=AlignerResultType.EVALUATED_PREDICTIONS
            )

            # Only the prediction with best score is kept (=BEST PREDICTION) -idea? also make a list of all predictions?
            for key, results in dict_predictions_evaluated.items():
                # call the parcels, evaluate through brdr and merge them
                rd = next(iter(results))
                dict_prediction_rd[key] = rd
                state = results[rd]["properties"]["brdr_evaluation"].value
                dict_prediction_state[key] = state
                score = results[rd]["properties"]["brdr_prediction_score"]
                dict_prediction_score[key] = score
                print(
                    f"brdr correction: '{state}' with relevant distance {rd}  & score {str(score)} for ID {key}"
                )
                geom = results[rd]["result"]
                dict_prediction_wkt[key] = geom_to_wkt(geom)
                brdr_parcels, parcel_100perc_list = get_parcel_lists(aligner, geom)
                brdr_parcel_count = len(brdr_parcels)
                dict_prediction[key] = brdr_parcel_count
                # print (brdr_parcel_count)

            fp_parcels = parcel_buffer_counts[1]
            # En welke zijn dit dan? Deze mee oplijsten
            if not (score is None or score == -1):
                # vip-buffer minus brdr
                fp_parcels = list(
                    set(buffer_data[0.2]["parcels"])
                    - set(brdr_parcels)
                    - set(coverage_data[100]["parcels"])
                )
                # calculated_false_positives = parcel_buffer_counts[1] - brdr_parcel_count
                s1, s2, s3 = (
                    set(buffer_data[2]["parcels"]),
                    set(coverage_data[50]["parcels"]),
                    set(brdr_parcels),
                )
                all_ids = s1 | s2 | s3
                common_ids = s1 & s2 & s3
            else:
                # vip-buffer minus buffer_2
                fp_parcels = list(
                    set(buffer_data[0.2]["parcels"])
                    - set(buffer_data[5]["parcels"])
                    - set(coverage_data[100]["parcels"])
                )
                # calculated_false_positives = parcel_buffer_counts[1] -parcel_buffer_counts[-1]
                s1, s2 = set(buffer_data[2]["parcels"]), set(
                    coverage_data[50]["parcels"]
                )
                all_ids = s1 | s2
                common_ids = s1 & s2

            doubt_parcels = list(all_ids - common_ids)

            # Do we want to add the extra doubtable parcels to the 'estimated false positive parcels'?
            # Currently decided to not do that -> line below commented
            # fp_parcels = list(set(fp_parcels) | set(doubt_parcels))
            fp_parcel_count = len(fp_parcels)
            fp_geometry = safe_unary_union(
                [aligner.reference_data.features[p].geometry for p in fp_parcels]
            ).wkt
            dict_fp[key] = fp_parcel_count
            dict_fp_parcels[key] = fp_parcels
            dict_fp_wkt[key] = fp_geometry
            # search for ids where the different methods give another result

            dict_doubt_parcels[key] = doubt_parcels
            doubt_geometry = safe_unary_union(
                [aligner.reference_data.features[p].geometry for p in doubt_parcels]
            ).wkt
            dict_doubt_wkt[key] = doubt_geometry

            metrics["success_ids"].append(key)

        except Exception as e:
            if key is None:
                raise KeyError(f"Error - No key found")
            metrics["failed_ids"].append(key)
            print(f"Error processing feature {key}: {e}")
            dict_prediction_wkt[key] = GeometryCollection().wkt
            dict_fp_wkt[key] = GeometryCollection().wkt
            dict_doubt_wkt[key] = GeometryCollection().wkt

    # Combine in dataframe
    df = pd.DataFrame(
        {
            "coverage_analysis": dict_coverage_list,
            "coverage_range": dict_coverage_range,
            "coverage_0": dict_coverage_0,
            "coverage_50": dict_coverage_50,
            "buffer_analysis": dict_buffer_list,
            "buffer_vip": dict_buffer_vip,
            "buffer_2": dict_buffer_2,
            "buffer_5": dict_buffer_5,
            FP_ESIMATION_COLUMN_NAME: dict_fp,
            "brdr_prediction_5m": dict_prediction,
            "brdr_rd": dict_prediction_rd,
            "brdr_state": dict_prediction_state,
            "brdr_score": dict_prediction_score,
            "fp_parcels": dict_fp_parcels,
            "doubt_parcels": dict_doubt_parcels,
            DOUBT_WKT: dict_doubt_wkt,
            FALSE_POSITIVE_WKT: dict_fp_wkt,
            BRDR_WKT: dict_prediction_wkt,
            ORIGINAL_WKT: dict_original_wkt,
        }
    )

    # Set index to column 'id'
    df.reset_index(inplace=True)
    df.rename(columns={"index": id_name}, inplace=True)
    now = datetime.now()
    metrics["endtime"] = now
    metrics["total_time"] = now - starttime

    return df, metrics


def get_folder_path(analysis_name):
    # Generate timestamp
    foldername = analysis_name + "_" + datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = Path(foldername)
    output_dir.mkdir(parents=True, exist_ok=True)

    return output_dir


def export_wkt_columns_to_geojson(df, wkt_columns, path, crs=DEFAULT_CRS):
    """
    Exports multiple WKT columns from a DataFrame to separate GeoJSON files,
    ensuring that other WKT columns are excluded from the properties.

    Parameters
    ----------
    df : pd.DataFrame
        The source DataFrame containing geometry and property columns.
    wkt_columns : list of str
        The names of the columns containing WKT strings.
    path : str or Path
        The directory where the GeoJSON files will be saved.
    crs : str, optional
        The coordinate reference system for the GeoJSON.
    """
    output_path = Path(path)
    output_path.mkdir(parents=True, exist_ok=True)

    for col in wkt_columns:
        if col not in df.columns:
            print(f"Warning: Column '{col}' not found in DataFrame. Skipping.")
            continue

        # 1. Identify property columns: Exclude ALL columns listed in wkt_columns
        # This ensures that no WKT strings remain in the attributes
        property_cols = [c for c in df.columns if c not in wkt_columns]

        # 2. Select properties + the current WKT column
        temp_df = df[property_cols + [col]].copy()

        # 3. Convert WKT string to actual geometry
        temp_df["geometry"] = temp_df[col].apply(wkt.loads)

        # 4. Initialize GeoDataFrame
        gdf = gpd.GeoDataFrame(temp_df, geometry="geometry", crs=crs)

        # 5. Drop the redundant string version of the current column
        gdf = gdf.drop(columns=[col])

        # 6. Export to file
        file_name = f"{col}.geojson"
        final_destination = output_path / file_name

        gdf.to_file(final_destination, driver="GeoJSON")
        print(
            f"Successfully exported {col} (properties cleaned) to {final_destination}"
        )


def export_stats(df, column_name, tolerance, path, filename=STATS_TXT):
    """
    Calculates error statistics for a specific column and exports them to a text file.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame containing the data to analyze.
    column_name : str
        The name of the column used for statistics (e.g., 'error_distance').
    tolerance : float
        The threshold value to calculate the percentage of values within limits.
    path : str or Path
        The file path where the stats.txt will be saved.

    Returns
    -------
    None
    """
    if column_name not in df.columns:
        raise ValueError(f"Column '{column_name}' not found in DataFrame.")

        # 1. Perform calculations
    data = df[column_name].dropna()
    if data.empty:
        print(f"Warning: Column '{column_name}' is empty after dropping NaNs.")
        return

    mean_error = data.mean()
    median_error = data.median()
    max_error = data.max()

    # RMS Calculation
    rms_error = np.sqrt(np.mean(data**2))

    # Percentiles calculation (P25, P50, P75, P90, P95)
    percentiles = data.quantile([0.25, 0.5, 0.75, 0.9, 0.95]).to_dict()

    # Percentage within tolerance
    within_tol_count = (data <= tolerance).sum()
    within_tol_percent = (within_tol_count / len(data)) * 100

    # 2. Prepare stats dictionary
    stats = {
        "Mean error": mean_error,
        "Median error": median_error,
        "Max error": max_error,
        "RMS error": rms_error,
        f"% beneath threshold ({tolerance})": within_tol_percent,
        "---": None,  # Visual separator
        "P25 (1st Quartile)": percentiles[0.25],
        "P50 (Median)": percentiles[0.5],
        "P75 (3rd Quartile)": percentiles[0.75],
        "P90": percentiles[0.9],
        "P95": percentiles[0.95],
    }
    # 3. Write to file
    file_path = Path(path / filename)

    with open(file_path, "w", encoding="utf-8") as f:
        f.write(f"Stats for column: {column_name}\n")
        f.write("-" * 40 + "\n")
        for k, v in stats.items():
            if v is None:
                # Seperator
                f.write("-" * 45 + "\n")
            else:
                # Format numbers
                f.write(f"{k:<30}: {v:.3f}\n")

    print(f"Statistics for '{column_name}' written to: {file_path}")


def export_metrics(metrics, path, filename=METRICS_TXT):
    """Internal helper to write the counter dictionary to a text file."""
    path = Path(path)
    dataset_date = str(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    file_path = Path(path / filename)

    with open(file_path, "a", encoding="utf-8") as f:
        f.write("\nDetailed Processing Report\n")
        f.write("=" * 45 + "\n")
        f.write(f"\nStart time: {metrics['starttime']}\n")
        f.write(f"\nEnd time: {metrics['endtime']}\n")
        f.write(f"\nTotal time: {metrics['total_time']}\n")
        f.write(f"\nArea limit: {metrics['area_limit']}\n")
        f.write(f"\nNameDate: {dataset_date}\n")
        f.write(f"\nArea limit: {metrics['area_limit']}\n")
        f.write("=" * 45 + "\n")
        f.write(f"{'Total features loaded':<30}: {metrics['total_count']}\n")
        f.write(f"{'Successfully processed':<30}: {len(metrics['success_ids'])}\n")
        f.write(f"{'Skipped by area limit' :<30}: {len(metrics['skipped_ids'])}\n")
        f.write(f"{'Failed (Errors)':<30}: {len(metrics['failed_ids'])}\n")

        if metrics["skipped_ids"]:
            f.write(
                "\nSkipped IDs:\n" + ", ".join(map(str, metrics["skipped_ids"])) + "\n"
            )

        if metrics["failed_ids"]:
            f.write(
                "\nFailed IDs:\n" + ", ".join(map(str, metrics["failed_ids"])) + "\n"
            )
        f.write("-" * 45 + "\n")


def export_analysis_results(
    path,
    df=None,
    metrics=None,
    column_name=FP_ESIMATION_COLUMN_NAME,
    tolerance=2,
    wkt_columns=[ORIGINAL_WKT, BRDR_WKT, FALSE_POSITIVE_WKT, DOUBT_WKT],
):
    if df is None:
        df = pd.read_csv(path / FP_ANALYSIS_CSV_NAME)
    if df is None or df.empty:
        print("Warning: No data to export.")
        return
    # Export while dropping wkt columns
    df.drop(columns=wkt_columns).to_csv(path / FP_ANALYSIS_CSV_NAME, index=False)
    if metrics:
        export_metrics(metrics=metrics, path=path, filename=METRICS_TXT)
    export_stats(df=df, column_name=column_name, tolerance=tolerance, path=path)
    export_wkt_columns_to_geojson(df=df, wkt_columns=wkt_columns, path=path)
    export_boxplot(df=df, path=path, column_name=column_name)
    export_histogram(df=df, path=path, column_name=column_name)
