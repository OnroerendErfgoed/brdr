import logging
from math import ceil

import geopandas as gpd
import matplotlib.pyplot as plt

from brdr.typings import ProcessResult


def _make_map(ax, processresult, thematic_dict, reference_dict):
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
        dicts = _processresult_to_dicts(processresult)
        results = dicts[0]
        results_diff_pos = dicts[1]
        results_diff_neg = dicts[2]
        if ax is None:
            ax = plt.subplot(1, 1, 1)
        # ax_result =
        gpd.GeoSeries(list(results.values())).plot(
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
        # ax_diff_pos = (
        gpd.GeoSeries(list(results_diff_pos.values())).plot(
            ax=ax,
            color="none",
            edgecolor="green",
            hatch="+",
            linewidth=0.0,
            linestyle="dashdot",
            label="diff_plus",
            zorder=4,
        )
        # ax_diff_neg =
        gpd.GeoSeries(list(results_diff_neg.values())).plot(
            ax=ax,
            color="none",
            edgecolor="red",
            hatch="+",
            linewidth=0.0,
            linestyle="dashdot",
            label="diff_min",
            zorder=5,
        )
        # save the extent of original, resulting and difference - geometries
        axis_extent = list(ax_thematic_dict.viewLim.intervalx) + list(
            ax_thematic_dict.viewLim.intervaly
        )
        # ax_reference_dict =
        gpd.GeoSeries(list(reference_dict.values())).plot(
            ax=ax,
            color="#FFF8C9",
            edgecolor="black",
            linewidth=2.0,
            label="reference",
            zorder=1,
        )
        # zoom map to saved extent
        ax.axis(axis_extent)
    except Exception:  # noqa
        logging.error("make_map: Error while making map")
    return ax


def show_map(
    dict_results: dict[str, dict[float, ProcessResult]],
    dict_thematic,
    dict_reference,
):
    """
    Show results on a map
    """
    dict_results_by_distance = {}
    for theme_id, dist_result in dict_results.items():
        for rel_dist, processresults in dist_result.items():
            dict_results_by_distance[rel_dist] = {}
            dict_results_by_distance[rel_dist][theme_id] = processresults

    len_series = len(dict_results_by_distance.keys())
    i = 0
    # Plot data in subplots
    len_series_half = ceil(len_series / 2)  # calculate half of the length of the series
    for dist in dict_results_by_distance:
        ax = plt.subplot(len_series_half, 2, i + 1)
        ax = _make_map(
            ax,  # noqa
            dict_results_by_distance[dist],
            dict_thematic,
            dict_reference,
        )
        ax.set_title("Relevant distance (m):" + str(dist))
        i = i + 1
    # Adjust layout
    # plt.tight_layout()
    # Show figure
    plt.show()


def print_brdr_formula(dict_results, aligner):
    for theme_id, dist_results in dict_results.items():
        for rel_dist, processresults in dist_results.items():
            print(
                "--------Formula for ID  "
                + str(theme_id)
                + " with relevant distance "
                + str(rel_dist)
                + "--------------"
            )
            print(aligner.get_brdr_formula(processresults["result"]))
    return


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


def _processresult_to_dicts(processresult):
    """
    Transforms a dictionary with all ProcessResults to individual dictionaries of the
    results
    Args:
        processresult:

    Returns:

    """
    results = {}
    results_diff = {}
    results_diff_plus = {}
    results_diff_min = {}
    results_relevant_intersection = {}
    results_relevant_diff = {}
    for key in processresult:
        processresult = processresult[key]
        results[key] = processresult["result"]
        results_diff[key] = processresult["result_diff"]
        results_diff_plus[key] = processresult["result_diff_plus"]
        results_diff_min[key] = processresult["result_diff_min"]
        results_relevant_intersection[key] = processresult[
            "result_relevant_intersection"
        ]
        results_relevant_diff[key] = processresult["result_relevant_diff"]

    return (
        results,
        results_diff,
        results_diff_plus,
        results_diff_min,
        results_relevant_intersection,
        results_relevant_diff,
    )
