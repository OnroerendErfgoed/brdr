from datetime import datetime
from math import inf
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from shapely import wkt

from brdr.aligner import Aligner
from brdr.be.grb.enums import GRBType
from brdr.be.grb.loader import GRBActualLoader
from brdr.configs import AlignerConfig
from brdr.constants import RELEVANT_DISTANCE_DECIMALS, DEFAULT_CRS
from brdr.enums import AlignerResultType
from brdr.geometry_utils import buffer_neg, geom_to_wkt, safe_unary_union
from brdr.loader import DictLoader
from brdr.utils import geojson_geometry_to_shapely

FALSE_POSITIVE_WKT = "fp_wkt"
DOUBT_WKT = "doubt_wkt"
BRDR_WKT="brdr_wkt"

FP_ANALYSIS_CSV_NAME = "fp_analysis.csv"
FP_HISTOGRAM_NAME = "fp_histogram.png"
FP_BOXPLOT_NAME = "fp_boxplot.png"
FP_ESIMATION_COLUMN_NAME= "fp_estimation"

def get_parcel_lists(aligner,geom):
    parcel_100perc_list = []
    observation = aligner.compare_to_reference(geom)
    reference_features_list=[]
    if "reference_features" in observation:
        reference_features_list=observation["reference_features"].keys()
        for k,v in observation["reference_features"].items():
            if v["full"]:
                parcel_100perc_list.append(k)
    return list(reference_features_list),parcel_100perc_list


def get_buffer_data(aligner, geom, buffers = [0.1, 0.2, 0.5, 1, 2, 3, 4, 5]):
    buffer_dict={}
    for buffer in buffers:
        buffer_dict[buffer] = {}
        geom_buffered = buffer_neg(geom,buffer)
        parcel_list,parcel_100perc_list =get_parcel_lists(aligner,geom_buffered)
        buffer_dict[buffer] ["geometry"]=geom_buffered
        buffer_dict[buffer] ["parcels"]=parcel_list
        buffer_dict[buffer] ["parcels_full"]=parcel_100perc_list

    return buffer_dict

def get_coverage_data(aligner, geom, percentages = [0, 1, 5, 10, 50, 90, 95, 99, 100]):

    coverage_dict = {item: {"parcels":[],"count":0} for item in percentages}

    observation = aligner.compare_to_reference(geom)
    if "reference_features" in observation:
        for k,v in observation["reference_features"].items():
            area_percentage =v["percentage"]
            for p in percentages:
                if area_percentage>=p:
                    coverage_dict[p]["parcels"].append(k)
                    coverage_dict[p]["count"]+=1

    return coverage_dict

def get_false_positive_grb_parcels_dataframe(data, id_name, area_limit=inf, processor=None):
    rds = [
        round(k, RELEVANT_DISTANCE_DECIMALS)
        for k in np.arange(0, 510, 10, dtype=int) / 100
    ]
    dict_coverage_list ={}
    dict_coverage_range ={}
    dict_coverage_0 ={}
    dict_coverage_50 ={}
    dict_buffer_list ={}
    dict_buffer_vip = {}
    dict_buffer_2 = {}
    dict_buffer_5 = {}
    dict_prediction_rd = {}
    dict_prediction = {}
    dict_prediction_wkt = {}
    dict_prediction_state={}
    dict_prediction_score={}
    dict_fp ={}
    dict_fp_parcels = {}
    dict_fp_wkt = {}
    dict_doubt_parcels = {}
    dict_doubt_wkt = {}
    features = []
    if data and "features" in data:
        features = data['features']
    print(f"Nr of  features: {len(features)}")
    for f in features:
        key = f["properties"][id_name]
        geom = geojson_geometry_to_shapely(f["geometry"])
        if geom.area > area_limit:
            print(f"Feature {key} skipped based on area {geom.area} m²")
            continue
        print(key)
        dict_theme = {key: geom}

        aligner_config = AlignerConfig
        aligner = Aligner(crs="EPSG:31370", processor=processor, config=aligner_config)
        loader = DictLoader(dict_theme)
        aligner.load_thematic_data(loader)
        loader = GRBActualLoader(
            grb_type=GRBType.ADP, partition=1000, aligner=aligner
        )
        aligner.load_reference_data(loader)


        # COVERAGE ANALYSIS
        coverage_data = get_coverage_data(aligner, geom, percentages=[0, 1, 5, 10, 50, 90, 95, 99, 100])
        parcel_coverage_counts = [coverage_data[b]["count"] for b in coverage_data.keys()]
        dict_coverage_list[key] = parcel_coverage_counts
        dict_coverage_range[key] = parcel_coverage_counts[0] - parcel_coverage_counts[-1]
        dict_coverage_0[key] = parcel_coverage_counts[0]
        dict_coverage_50[key] = parcel_coverage_counts[4]

        # BUFFER ANALYSIS
        buffer_data = get_buffer_data(aligner, geom, buffers=[0.1, 0.2, 0.5, 1, 2, 3, 4, 5])
        parcel_buffer_counts = [len(buffer_data[b]["parcels"]) for b in buffer_data.keys()]
        dict_buffer_list[key] = parcel_buffer_counts
        dict_buffer_vip[key] = parcel_buffer_counts[1]
        #parcels_vip = buffer_data.get(0.2)["parcels"]
        dict_buffer_2[key] = parcel_buffer_counts[4]
        dict_buffer_5[key] = parcel_buffer_counts[-1]

        # BRDR ANALYSIS
        # This is done for every feature seperately
        aligner_result = aligner.evaluate(
            relevant_distances=rds,
            max_predictions=1,
            multi_to_best_prediction=True
        )
        dict_predictions_evaluated = aligner_result.get_results(aligner=aligner,
                                                                result_type=AlignerResultType.EVALUATED_PREDICTIONS)

        # TODO? also make a list of all predictions?
        for key, results in dict_predictions_evaluated.items():
            # call the parcels, evaluate through brdr and merge them
            rd = next(iter(results))
            dict_prediction_rd[key] = rd
            state = results[rd]["properties"]["brdr_evaluation"].value
            dict_prediction_state[key] = state
            score = results[rd]["properties"]["brdr_prediction_score"]
            dict_prediction_score[key] = score
            print(f"brdr correction: '{state}' with relevant distance {rd}  & score {str(score)} for ID {key}")
            geom = results[rd]['result']
            dict_prediction_wkt[key] = geom_to_wkt(geom)
            brdr_parcels, parcel_100perc_list = get_parcel_lists(aligner, geom)
            brdr_parcel_count = len(brdr_parcels)
            dict_prediction[key] = brdr_parcel_count
            # print (brdr_parcel_count)

        fp_parcels = parcel_buffer_counts[1]
        # En welke zijn dit dan? Deze mee oplijsten
        if not (score is None or score == -1):
            # vip-buffer minus brdr
            fp_parcels = list(set(buffer_data[0.2]["parcels"]) - set(brdr_parcels) - set(coverage_data[100]["parcels"]))
            # calculated_false_positives = parcel_buffer_counts[1] - brdr_parcel_count
            s1, s2, s3 = set(buffer_data[2]["parcels"]), set(coverage_data[50]["parcels"]), set(brdr_parcels)
            all_ids = s1 | s2 | s3
            common_ids = s1 & s2 & s3
        else:
            # vip-buffer minus buffer_2
            fp_parcels = list(
                set(buffer_data[0.2]["parcels"]) - set(buffer_data[5]["parcels"]) - set(coverage_data[100]["parcels"]))
            # calculated_false_positives = parcel_buffer_counts[1] -parcel_buffer_counts[-1]
            s1, s2 = set(buffer_data[2]["parcels"]), set(coverage_data[50]["parcels"])
            all_ids = s1 | s2
            common_ids = s1 & s2

        doubt_parcels = list(all_ids - common_ids)
        # fp_parcels = list(set(fp_parcels) | set(doubt_parcels)) #TODO or not?
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
        }
    )

    # Set index to column 'id'
    df.reset_index(inplace=True)
    df.rename(columns={"index": id_name}, inplace=True)
    return df


def get_folder_path(analysis_name):
    # Generate timestamp
    foldername = analysis_name + "_" + datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = Path(foldername)
    output_dir.mkdir(parents=True, exist_ok=True)

    return output_dir

def export_wkt_columns_to_geojson(df, wkt_columns, path):
    """
    Exports multiple WKT columns from a DataFrame to separate GeoJSON files.

    Parameters
    ----------
    df : pd.DataFrame
        The source DataFrame containing geometry and property columns.
    wkt_columns : list of str
        The names of the columns containing WKT strings.
    output_dir : str or Path
        The directory where the GeoJSON files will be saved.
    """
    output_path = Path(path)
    output_path.mkdir(parents=True, exist_ok=True)

    for col in wkt_columns:
        if col not in df.columns:
            print(f"Warning: Column '{col}' not found in DataFrame. Skipping.")
            continue

        # 1. Create a copy to avoid modifying the original DF
        # We only take the geometry column and the other non-WKT columns as properties
        other_cols = [c for c in df.columns if c not in wkt_columns]
        temp_df = df[other_cols + [col]].copy()

        # 2. Convert WKT strings to actual geometry objects
        temp_df["geometry"] = temp_df[col].apply(wkt.loads)

        # 3. Create GeoDataFrame
        gdf = gpd.GeoDataFrame(temp_df, geometry="geometry", crs=DEFAULT_CRS)

        # 4. Remove the original WKT string column to avoid redundancy in properties
        gdf = gdf.drop(columns=[col])

        # 5. Define filename and export
        file_name = f"{col}.geojson"
        final_destination = output_path / file_name

        gdf.to_file(final_destination, driver="GeoJSON")
        print(f"Successfully exported {col} to {final_destination}")


def export_stats(df, column_name, tolerance, path,filename="stats.txt"):
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
    data = df[column_name].dropna()  # Exclude missing values for accuracy

    mean_error = data.mean()
    max_error = data.max()

    # RMS Calculation: Square root of the mean of the squares
    rms_error = np.sqrt(np.mean(data**2))

    # Percentage within tolerance
    within_tol_count = (data <= tolerance).sum()
    within_tol_percent = (within_tol_count / len(data)) * 100 if len(data) > 0 else 0

    # 2. Prepare stats dictionary
    stats = {
        "Mean error": mean_error,
        "Max error": max_error,
        "RMS error": rms_error,
        f"% beneath threshold ({tolerance})": within_tol_percent,
    }

    # 3. Write to file
    file_path = Path(path /filename)

    with open(file_path, "w", encoding="utf-8") as f:
        f.write(f"Stats for column: {column_name}\n")
        f.write("-" * 40 + "\n")
        for k, v in stats.items():
            # Format numbers to 3 decimal places
            f.write(f"{k:<30}: {v:.3f}\n")

    print(f"Statistics for '{column_name}' written to: {file_path}")

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def export_histogram(path, df, column_name, filename="histogram.png"):
    """
    Creates and saves a high-quality histogram from a DataFrame column.

    Parameters
    ----------
    path : str or Path
        Directory where the plot will be saved.
    df : pd.DataFrame
        The DataFrame containing the data.
    column_name : str
        The name of the column to plot.
    filename : str, optional
        The name of the resulting image file. Default is "histogram.png".

    Returns
    -------
    None
    """
    if df is None or df.empty:
        print("Warning: No data available to plot histogram.")
        return

    values = df[column_name].dropna().values
    if len(values) == 0:
        print("Warning: Column is empty. Skipping histogram.")
        return

    # Set modern style
    plt.style.use("seaborn-v0_8-muted")
    fig, ax = plt.subplots(figsize=(10, 6))

    # Define bins (integer steps)
    bins = np.arange(int(min(values)), int(max(values)) + 2, 1)

    # Plot histogram
    ax.hist(
        values,
        bins=bins,
        color="royalblue",
        edgecolor="white",
        linewidth=1.2,
        alpha=0.85,
        label="Object Count",
    )

    # X-axis configuration
    ticks = np.arange(int(min(values)), int(max(values)) + 1, step=2)
    ax.set_xticks(ticks)

    # Styling
    ax.yaxis.grid(True, linestyle="--", alpha=0.7)
    ax.set_axisbelow(True)
    ax.set_title(f"Histogram - {column_name}", fontsize=14, pad=15, fontweight="bold")
    ax.set_xlabel("Value", fontsize=12)
    ax.set_ylabel("Count", fontsize=12)

    # Remove top and right spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.legend()
    plt.tight_layout()

    save_path = Path(path) / filename
    plt.savefig(save_path, dpi=300)
    plt.close(fig)  # Close to free up memory
    print(f"Histogram saved to: {save_path}")


def export_boxplot(path, df, column_name, filename="boxplot.png"):
    """
    Creates and saves a horizontal boxplot from a DataFrame column.

    Parameters
    ----------
    path : str or Path
        Directory where the plot will be saved.
    df : pd.DataFrame
        The DataFrame containing the data.
    column_name : str
        The name of the column to plot.
    filename : str, optional
        The name of the resulting image file. Default is "boxplot.png".

    Returns
    -------
    None
    """
    if df is None or df.empty:
        print("Warning: No data available to plot boxplot.")
        return

    values = df[column_name].dropna().values
    if len(values) == 0:
        print("Warning: Column is empty. Skipping boxplot.")
        return

    plt.style.use("seaborn-v0_8-muted")
    fig, ax = plt.subplots(figsize=(10, 4))

    # Plot boxplot
    result = ax.boxplot(
        values,
        vert=False,
        patch_artist=True,
        notch=False,
        showmeans=True,
        meanprops={
            "marker": "^",
            "markerfacecolor": "green",
            "markeredgecolor": "green",
        },
        medianprops={"color": "black", "linewidth": 2},
        flierprops={
            "markerfacecolor": "red",
            "marker": "o",
            "markersize": 5,
            "alpha": 0.5,
        },
    )

    # Coloring the box
    for patch in result["boxes"]:
        patch.set_facecolor("royalblue")
        patch.set_alpha(0.7)

    # X-axis configuration
    ticks = np.arange(int(min(values)), int(max(values)) + 1, step=2)
    ax.set_xticks(ticks)

    # Styling
    ax.set_title(f"Boxplot - {column_name}", fontsize=14, fontweight="bold", pad=15)
    ax.set_xlabel("Value", fontsize=12)
    ax.set_yticks([])  # Hide Y-axis labels for a single box

    ax.grid(axis="x", linestyle="--", alpha=0.6)
    ax.set_axisbelow(True)

    # Clean up spines
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)

    plt.tight_layout()

    save_path = Path(path) / filename
    plt.savefig(save_path, dpi=300)
    plt.close(fig)
    print(f"Boxplot saved to: {save_path}")

def export_analysis_results(path,df=None, column_name=FP_ESIMATION_COLUMN_NAME, tolerance=2,wkt_columns=[BRDR_WKT, FALSE_POSITIVE_WKT, DOUBT_WKT]):
    if df is None:
        df = pd.read_csv(path / FP_ANALYSIS_CSV_NAME)
    if df is None or df.empty:
        print("Warning: No data available to plot histogram.")
        return
    export_stats(df=df, column_name=column_name, tolerance=tolerance, path=path)
    export_wkt_columns_to_geojson(df=df, wkt_columns=wkt_columns, path=path)
    export_boxplot(df=df, path=path, column_name=column_name)
    export_histogram(df=df, path=path, column_name=column_name)