import logging
from datetime import datetime
from math import ceil
from pathlib import Path
from pprint import pprint

import geopandas as gpd
import numpy as np
from shapely.geometry import Point, LineString

from brdr.constants import DEFAULT_CRS
from brdr.typings import ProcessResult

try:
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter
    from matplotlib.patches import Patch
    from PIL import Image


except ImportError:
    raise ImportError(
        "The visualisation-functions require 'matplotlib' & 'pillow'. "
        "Install these by: pip install brdr[viz]"
    )


def _make_map(ax, processresult, thematic_dict, reference_dict):
    """
    Render a map zoomed strictly to the bounds of thematic geometries.
    """
    try:
        dicts = _processresult_to_dicts(processresult)
        results = dicts[0]
        results_diff_pos = dicts[2]
        results_diff_neg = dicts[3]

        if ax is None:
            ax = plt.subplot(1, 1, 1)

        # We maken GeoSeries van de data
        gs_thematic = gpd.GeoSeries(list(thematic_dict.values()))
        gs_reference = gpd.GeoSeries(list(reference_dict.values()))
        gs_results = gpd.GeoSeries(list(results.values()))

        # --- FOCUS LOGICA ---
        # Bereken de exacte bounding box van de thematische data: [minx, miny, maxx, maxy]
        bounds = gs_thematic.total_bounds
        padding = 2.0  # Voeg 2 meter (of graden, afhankelijk van CRS) marge toe

        # Stel de limieten in VÓÓR of NA het plotten, dit overschrijft de automatische schaal
        ax.set_xlim(bounds[0] - padding, bounds[2] + padding)
        ax.set_ylim(bounds[1] - padding, bounds[3] + padding)

        # 1. Plot Reference Data (wordt nu 'afgesneden' door bovenstaande limieten)
        gs_reference.plot(
            ax=ax, color="#FFF8C9", edgecolor="black", linewidth=2.0, zorder=1
        )

        # 2. Plot Resulting Geometry
        gs_results.plot(
            ax=ax, alpha=0.8, color="none", edgecolor="green", linewidth=5.0, zorder=2
        )

        # 3. Plot Original (Thematic) Data
        gs_thematic.plot(
            ax=ax,
            alpha=0.6,
            color="none",
            edgecolor="#0000FF",
            linewidth=3.0,
            linestyle="dashdot",
            zorder=3,
        )

        # 4. Plot Differences
        gpd.GeoSeries(list(results_diff_pos.values())).plot(
            ax=ax, color="none", edgecolor="green", hatch="+", zorder=4
        )
        gpd.GeoSeries(list(results_diff_neg.values())).plot(
            ax=ax, color="none", edgecolor="red", hatch="+", zorder=5
        )

        # Legenda
        legend_patch_res = Patch(facecolor="none", edgecolor="green", label="Result")
        legend_patch_ori = Patch(
            facecolor="none", edgecolor="blue", label="Original", linestyle="dashdot"
        )
        ax.legend(handles=[legend_patch_res, legend_patch_ori])

    except Exception as e:
        logging.error(f"make_map: Error while making map: {e}")

    return ax


def show_map(
    aligner_results: dict[any, dict[float, ProcessResult]],
    dict_thematic,
    dict_reference,
):
    """
    Display a grid of maps showing alignment results across different distances.

    This function creates a multi-panel figure where each subplot represents
    the geographic state of the thematic data after being aligned using a
    specific 'relevant distance' parameter.

    Parameters
    ----------
    aligner_results : dict
        A nested dictionary structured as `{theme_id: {rel_dist: ProcessResult}}`.
    dict_thematic : dict
        Dictionary containing the original thematic geometries.
    dict_reference : dict
        Dictionary containing the reference geometries used for alignment.

    Returns
    -------
    None
        The function generates and displays a Matplotlib figure.

    Notes
    -----
    The function automatically calculates the required number of rows for
    a two-column grid. It relies on the internal `_make_map` helper
    function to render the individual spatial layers.
    """
    # Re-index to group by distance instead of ID
    dict_results_by_distance = {}
    for theme_id, dist_results in aligner_results.items():
        for rel_dist, processresults in dist_results.items():
            if rel_dist not in dict_results_by_distance:
                dict_results_by_distance[rel_dist] = {}
            dict_results_by_distance[rel_dist][theme_id] = processresults

    distances = sorted(dict_results_by_distance.keys())
    num_plots = len(distances)
    num_cols = 2
    num_rows = ceil(num_plots / num_cols)

    # Create figure with dynamic height
    # plt.figure(figsize=(12, 4 * num_rows))
    fig, ax = plt.subplots(figsize=(12, 4 * num_rows))

    for i, dist in enumerate(distances):
        ax = plt.subplot(num_rows, num_cols, i + 1)

        # Draw the map using the helper
        _make_map(
            ax,
            dict_results_by_distance[dist],
            dict_thematic,
            dict_reference,
        )

        ax.set_title(f"Relevant distance: {dist} m")

    plt.show(block=False)


def animated_map(
    aligner_results,
    dict_thematic,
    dict_reference,
    xlim,
    ylim,
    interval,
    filename,
):
    """
    Create a synchronized animation of map changes and area statistics.

    This function generates a side-by-side animation: the left panel shows the
    geometric alignment on a map (zoomed to thematic data), while the right
    panel plots the percentage of area change as a function of the
    'relevant distance' parameter.

    Parameters
    ----------
    aligner_results : dict
        A nested dictionary structured as `{theme_id: {rel_dist: ProcessResult}}`.
    dict_thematic : dict
        Dictionary containing the original thematic geometries.
    dict_reference : dict
        Dictionary containing the reference geometries used for alignment.
    xlim : float
        The maximum value for the X-axis (relevant distance) in the trend plot.
    ylim : float
        The maximum value for the Y-axis (% area change) in the trend plot.
    interval : int
        Delay between frames in milliseconds.
    filename : str
        The output path/name for the saved animation (e.g., 'analysis.gif').

    Returns
    -------
    None
        The function saves the animation to the specified `filename` and
        displays the plot.

    Notes
    -----
    To prevent "stray lines" in the trend plot during animation loops, the
    data tracking lists are cleared whenever the animation restarts at the
    first distance frame.
    """
    # Re-index results by distance for sequential animation
    dict_results_by_distance = {}
    for theme_id, dist_result in aligner_results.items():
        for rel_dist, processresults in dist_result.items():
            if rel_dist not in dict_results_by_distance:
                dict_results_by_distance[rel_dist] = {}
            dict_results_by_distance[rel_dist][theme_id] = processresults

    distances = sorted(list(dict_results_by_distance.keys()))

    # Figure setup: Map on left, Plot on right
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Initialize trend plot data tracking
    x_vals = []
    y_vals = []
    (plot_line,) = ax2.plot([], [], color="blue", lw=2, label="% area change")
    marker_line = ax2.axvline(
        x=0, color="red", linestyle="--", label="Current Relevant Distance"
    )

    ax2.set_xlim(0, xlim)
    ax2.set_ylim(0, ylim)
    ax2.set_title("Area change by relevant distance")
    ax2.set_xlabel("Relevant Distance")
    ax2.set_ylabel("Difference (% area)")
    ax2.legend(loc="upper left")
    ax2.grid(True, alpha=0.3)

    def update(frame_dist):
        """Internal update function for Matplotlib FuncAnimation."""

        # FIX: Clear lists if we are back at the first frame to avoid stray lines
        if frame_dist == distances[0]:
            x_vals.clear()
            y_vals.clear()

        # Update Map Subplot
        ax1.clear()
        _make_map(
            ax1,
            dict_results_by_distance[frame_dist],
            dict_thematic,
            dict_reference,
        )

        # Calculate stats for the current frame
        # We assume the first theme_id is the primary subject for the trend
        theme_id = list(dict_results_by_distance[frame_dist].keys())[0]
        entry = dict_results_by_distance[frame_dist][theme_id]

        area_diff = entry["result_diff"].area
        area_total = entry["result"].area
        percentage = 100 * area_diff / area_total if area_total != 0 else 0

        # Append new data points
        x_vals.append(frame_dist)
        y_vals.append(percentage)

        # Update trend plot elements
        plot_line.set_data(x_vals, y_vals)
        marker_line.set_xdata([frame_dist])

        return ax1, plot_line, marker_line

    # Create the animation object
    ani = FuncAnimation(fig, update, frames=distances, interval=interval, blit=False)

    # Save before showing (showing can sometimes clear the buffer)
    if filename:
        try:
            ani.save(filename, writer="pillow")
            print(f"Animation successfully saved to: {filename}")
        except Exception as e:
            logging.error(f"Failed to save animation: {e}")

    plt.tight_layout()
    plt.show(block=False)


def print_observation_of_aligner_results(aligner_results, aligner):
    """
    Print a detailed comparison analysis for alignment results.

    Iterates through a nested results dictionary and uses the Aligner's
    comparison logic to output a human-readable summary of how the
    processed geometry relates to the reference data for different
    parameter settings.

    Parameters
    ----------
    aligner_results : dict of dict
        A nested dictionary where the first key is the `theme_id` and
        the second key is the `relevant_distance` (rel_dist). The values
        are `ProcessResult` dictionaries.
    aligner : Aligner
        An instance of the Aligner class (or a subclass) that provides
        the `compare_to_reference` method for geometric analysis.

    Returns
    -------
    None
        The function prints the analysis directly to the standard output.

    Notes
    -----
    This function is primarily used for debugging and quality control
    to verify if the "Observation" (alignment logic) is producing the
    expected geometric relationship with the reference data.
    """
    for theme_id, dist_results in aligner_results.items():
        for rel_dist, processresults in dist_results.items():
            header = f"-------- Observation for ID {theme_id} with relevant distance {rel_dist} --------------"
            print(header)

            # Access the 'result' geometry and compare it to the reference
            comparison_report = aligner.compare_to_reference(processresults["result"])
            pprint(comparison_report)

    return


def plot_difference_by_relevant_distance(
    dict_differences,
    xlabel="relevant distance",
    ylabel="difference",
    title="Relevant distance vs difference",
    save_path=None,
):
    """
    Plot geometric differences across multiple features as a line graph.

    This function iterates through a nested dictionary where each key represents
    a feature (or scenario) and its value is a dictionary mapping X-values
    (typically 'relevant distance') to Y-values (geometric difference).

    Parameters
    ----------
    dict_differences : dict of dict
        A nested dictionary structure: `{label: {x_value: y_value}}`.
        Example: `{'feature_1': {0.1: 0.5, 0.2: 0.3}}`.
    xlabel : str, default "relevant distance"
        Label for the X-axis.
    ylabel : str, default "difference"
        Label for the Y-axis.
    title : str, default "Relevant distance vs difference"
        The main title of the plot.
    save_path : str or Path, optional
        If provided, the plot will be saved to this directory or file path.
        If a directory is provided, a timestamped filename is generated.

    Returns
    -------
    None
        The function displays the plot using `plt.show()` and optionally
        saves the figure.

    Notes
    -----
    This function is particularly useful for sensitivity analysis—showing how
    adjusting the `relevant_distance` parameter impacts the final
    geometric outcome.
    """
    # plt.figure(figsize=(10, 6))
    fig, ax = plt.subplots(figsize=(10, 6))

    for key, diffs in dict_differences.items():
        # Ensure values are sorted by x-value to prevent 'zig-zag' lines
        sorted_keys = sorted(diffs.keys())
        x_values = sorted_keys
        y_values = [diffs[k] for k in sorted_keys]

        plt.plot(x_values, y_values, marker="o", linestyle="-", label=str(key))

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.grid(True, alpha=0.3)

    if save_path:
        path = Path(save_path)
        if path.is_dir() or not path.suffix:
            path.mkdir(parents=True, exist_ok=True)
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            path = path / f"diff_plot_{timestamp}.png"

        plt.savefig(path)
        print(f"Plot saved to: {path}")

    plt.show(block=False)


def _processresult_to_dicts(process_results):
    """
    Transform a nested dictionary of ProcessResults into categorized dictionaries.

    This helper function pivots data from a key-based structure where each key
    contains a full ProcessResult, into multiple dictionaries grouped by result
    type (e.g., all differences, all intersections).

    Parameters
    ----------
    process_results : dict
        A dictionary where keys are identifiers (e.g., feature IDs) and values
        are dictionaries containing specific keys: "result", "result_diff",
        "result_diff_plus", "result_diff_min", "result_relevant_intersection",
        and "result_relevant_diff".

    Returns
    -------
    tuple of dict
        A tuple containing six dictionaries in the following order:
        - results: The primary processed geometries.
        - results_diff: Total geometric differences.
        - results_diff_plus: Positive differences (additions).
        - results_diff_min: Negative differences (removals).
        - results_relevant_intersection: Intersections with relevant reference data.
        - results_relevant_diff: Differences within the relevant search area.

    Notes
    -----
    This is an internal helper function used to prepare data for bulk conversion
    into GeoDataFrames or for visualization purposes.
    """
    results = {}
    results_diff = {}
    results_diff_plus = {}
    results_diff_min = {}
    results_relevant_intersection = {}
    results_relevant_diff = {}

    for key, entry in process_results.items():
        results[key] = entry["result"]
        results_diff[key] = entry["result_diff"]
        results_diff_plus[key] = entry["result_diff_plus"]
        results_diff_min[key] = entry["result_diff_min"]
        results_relevant_intersection[key] = entry["result_relevant_intersection"]
        results_relevant_diff[key] = entry["result_relevant_diff"]

    return (
        results,
        results_diff,
        results_diff_plus,
        results_diff_min,
        results_relevant_intersection,
        results_relevant_diff,
    )


def export_to_geopackage(G, output_filename="graph_output.gpkg", crs=DEFAULT_CRS):
    """
    Export a NetworkX graph to a GeoPackage file with separate layers for nodes and edges.

    This function assumes that the node identifiers are tuples representing
    coordinates (x, y).
    All attributes attached to nodes and edges in the NetworkX graph are
    preserved as columns in the resulting attribute tables.

    Parameters
    ----------
    G : networkx.Graph
        The NetworkX graph object to be exported. Nodes should be (x, y) tuples.
    output_filename : str, optional
        The path and name of the output GeoPackage file, by default "graph_output.gpkg".
    crs : str, optional
        Coordinate system

    Returns
    -------
    None
        The function saves the file to disk and prints a confirmation message.

    Notes
    -----
    Example usage:
    export_to_qgis(my_graph, "network.gpkg")
    The output GeoPackage will contain two layers:
    1. 'nodes': Point geometries representing the graph vertices.
    2. 'edges': LineString geometries representing the connections.
    """

    # 1. Export Nodes
    node_data = []
    for node, data in G.nodes(data=True):
        # We assume 'node' is a coordinate tuple (x, y)
        node_data.append({"node_id": str(node), "geometry": Point(node), **data})
    gdf_nodes = gpd.GeoDataFrame(node_data, crs=crs)

    # 2. Export Edges
    edge_data = []
    for u, v, data in G.edges(data=True):
        edge_data.append(
            {"u": str(u), "v": str(v), "geometry": LineString([u, v]), **data}
        )
    gdf_edges = gpd.GeoDataFrame(edge_data, crs=crs)

    # 3. Save as layers in a single GeoPackage
    gdf_nodes.to_file(output_filename, layer="nodes", driver="GPKG")
    gdf_edges.to_file(output_filename, layer="edges", driver="GPKG")

    print(f"Exported graph:  You can now add this to GIS software.")


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
