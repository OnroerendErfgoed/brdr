import logging
import geopandas as gpd
import matplotlib.pyplot as plt
from math import ceil

import numpy as np


def _make_map(
    ax, results, results_diff_pos, results_diff_neg, thematic_dict, reference_dict
):
    """
    Fills an ax with a map:
     * reference_dict
     * theme_dict
     * resulting geometry
     * plus_differences
     * min_differences
    , so it can be used in matplotlib
    """
    try:
        if ax is None:
            ax = plt.subplot(1, 1, 1)
        ax_result = gpd.GeoSeries(list(results.values())).plot(
            ax=ax,
            alpha=0.5,
            color="none",
            hatch=" ",
            edgecolor="green",
            linewidth=7.0,
            label="result",
            zorder=2,
        )
        ax_thematic_dict = gpd.GeoSeries(list(thematic_dict.values())).plot(
            ax=ax,
            alpha=0.8,
            color="none",
            hatch="/",
            edgecolor="#0000FF",
            linewidth=3.0,
            linestyle="dashdot",
            label="theme",
            zorder=3,
        )
        ax_diff_pos = gpd.GeoSeries(list(results_diff_pos.values())).plot(
            ax=ax,
            color="none",
            edgecolor="green",
            hatch="+",
            linewidth=0.0,
            linestyle="dashdot",
            label="diff_plus",
            zorder=4,
        )
        ax_diff_neg = gpd.GeoSeries(list(results_diff_neg.values())).plot(
            ax=ax,
            color="none",
            edgecolor="red",
            hatch="+",
            linewidth=0.0,
            linestyle="dashdot",
            label="diff_min",
            zorder=5,
        )
        # save the extent of original, resuling and difference - geometries
        axis_extent = list(ax_thematic_dict.viewLim.intervalx) + list(
            ax_thematic_dict.viewLim.intervaly
        )
        ax_reference_dict = gpd.GeoSeries(list(reference_dict.values())).plot(
            ax=ax,
            color="#FFF8C9",
            edgecolor="black",
            linewidth=2.0,
            label="reference",
            zorder=1,
        )
        # zoom map to saved extent
        ax.axis(axis_extent)
    except:
        logging.error("make_map: Error while making map")
    return ax


def show_map(dict_results_by_distance, dict_thematic, dict_reference):
    """
    Show results on a map
    """
    len_series = len(dict_results_by_distance.keys())
    i = 0
    # Plot data in subplots
    len_series_half = ceil(len_series / 2)  # calculate half of the length of the series
    for dist in dict_results_by_distance:
        ax = plt.subplot(len_series_half, 2, i + 1)
        ax = _make_map(
            ax,
            dict_results_by_distance[dist][0],
            dict_results_by_distance[dist][2],
            dict_results_by_distance[dist][3],
            dict_thematic,
            dict_reference,
        )
        ax.set_title("Relevant distance (m):" + str(dist))
        i = i + 1
    # Adjust layout
    # plt.tight_layout()
    # Show figure
    plt.show()


def plot_series(
    series,
    dictionary,
    xlabel="relevant distance",
    ylabel="difference",
    title="Relevant distance vs difference",
):
    for key in dictionary:
        if len(dictionary[key]) == len(series):
            lst_diffs = list(dictionary[key].values())
            plt.plot(series, lst_diffs, label=str(key))
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.show()
    return
