import logging
from math import ceil

import geopandas as gpd
import matplotlib.pyplot as plt
from PIL import Image
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import Patch
from shapely.plotting import plot_polygon, plot_line, plot_points

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
        results_diff_pos = dicts[2]
        results_diff_neg = dicts[3]
        if ax is None:
            ax = plt.subplot(1, 1, 1)
        # ax_result =
        gpd.GeoSeries(list(results.values())).plot(
            ax=ax,
            alpha=0.8,
            color="none",
            hatch=" ",
            edgecolor="green",
            linewidth=5.0,
            label="Result",
            zorder=2,
        )

        ax_thematic_dict = gpd.GeoSeries(list(thematic_dict.values())).plot(
            ax=ax,
            alpha=0.6,
            color="none",
            hatch="/",
            edgecolor="#0000FF",
            linewidth=3.0,
            linestyle="dashdot",
            label="Original",
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
            label="PlusDifference",
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
            label="MinDifference",
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
            label="Reference",
            zorder=1,
        )
        # zoom map to saved extent
        ax.axis(axis_extent)
        legend_patch_res = Patch(
            facecolor="none", edgecolor="green", linewidth=1.0, label="Result"
        )
        legend_patch_ori = Patch(
            facecolor="none",
            edgecolor="blue",
            linestyle="dashdot",
            linewidth=1.0,
            label="Original",
        )
        # legend_patch_dif_min = Patch( edgecolor='red', linewidth=1.0, label='Min')
        # legend_patch_dif_plus = Patch( edgecolor='green', linewidth=1.0, label='Plus')
        ax.legend(
            handles=[
                legend_patch_res,
                legend_patch_ori,
                # legend_patch_dif_min,
                # legend_patch_dif_plus
            ]
            #    ,title="Layers"
        )
    except Exception:  # noqa
        logging.error("make_map: Error while making map")
    return ax


def show_map(
    dict_results: dict[any, dict[float, ProcessResult]],
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


def animated_map(
    dict_results: dict[any, dict[float, ProcessResult]],
    dict_thematic,
    dict_reference,
    xlim,
    ylim,
    interval,
    filename,
):
    """
    Show results on a map
    """
    dict_results_by_distance = {}
    for theme_id, dist_result in dict_results.items():
        for rel_dist, processresults in dist_result.items():
            dict_results_by_distance[rel_dist] = {}
            dict_results_by_distance[rel_dist][theme_id] = processresults

    distances = list(dict_results_by_distance.keys())
    len_distances = len(distances)

    # Voorbereiding van de figuur
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))
    _make_map(
        ax1,  # noqa
        dict_results_by_distance[distances[0]],
        dict_thematic,
        dict_reference,
    )

    # Rechter subplot: grafiek van oppervlakte over tijd
    x_vals = []
    y_vals = []
    (plot_line,) = ax2.plot([], [], color="blue", label="% area change")
    marker_line = ax2.axvline(
        x=0, color="red", linestyle="--", label="Relevant distance"
    )
    ax2.set_xlim(0, xlim)
    ax2.set_ylim(0, ylim)  # vaste Y-as schaal
    ax2.set_title("Area difference by relevant distance")
    ax2.set_xlabel("Actual relevant distance")
    ax2.set_ylabel("Difference (% area change)")
    ax2.legend()

    # Updatefunctie voor animatie
    def update(frame, extra_data):
        # Bereken nieuwe grootte
        print(str(frame))
        ax1 = fig.axes[0]
        ax1.clear()
        _make_map(
            ax1,  # noqa
            dict_results_by_distance[frame],
            dict_thematic,
            dict_reference,
        )
        key = list(dict_results_by_distance[frame].keys())[0]
        area_result_diff = dict_results_by_distance[frame][key]["result_diff"].area
        area_result = dict_results_by_distance[frame][key]["result"].area
        area = 100 * area_result_diff / area_result
        # Bereken oppervlakte
        x_vals.append(frame)
        y_vals.append(area)

        # Update grafiek
        plot_line.set_data(x_vals, y_vals)
        # f= int (frame)
        marker_line.set_xdata(frame)

        return plot_line, marker_line

    # Maak de animatie
    extra_data = None
    fps = 1000 / interval
    ani = FuncAnimation(
        fig, update, frames=distances, fargs=(extra_data,), interval=interval, blit=True
    )

    # Opslaan als animated GIF
    ani.save(filename, writer=PillowWriter(fps=fps))

    print("✅ GIF succesvol opgeslagen")


def print_brdr_observation(dict_results, aligner):
    for theme_id, dist_results in dict_results.items():
        for rel_dist, processresults in dist_results.items():
            print(
                "--------Formula for ID  "
                + str(theme_id)
                + " with relevant distance "
                + str(rel_dist)
                + "--------------"
            )
            print(aligner.compare_to_reference(processresults["result"]))
    return


def plot_dict_diffs(
    dict_diffs,
    xlabel="relevant distance",
    ylabel="difference",
    title="Relevant distance vs difference",
):
    for key, diffs in dict_diffs.items():
        x_values = list(diffs.keys())
        y_values = list(diffs.values())
        plt.plot(x_values, y_values, label=str(key))
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    #
    # from datetime import datetime
    #
    # now = datetime.now()  # current date and time
    # date_time = now.strftime("%m_%d_%Y_%H_%M_%S")
    # name =str(date_time) + ".png"
    # plt.savefig(name)
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


def show_geometry(geometry, ax=None, show=True):
    if ax is None:
        fig, ax = plt.subplots()
    geometry_type = geometry.geom_type
    if geometry_type in ["Polygon", "MultiPolygon"]:
        plot_polygon(geometry, ax=ax)
    elif geometry_type in ["LineString", "MultiLineString"]:
        plot_line(geometry, ax=ax)
    elif geometry_type in ["Point", "MultiPoint"]:
        plot_points(geometry, ax=ax)
    elif geometry_type in ["GeometryCollection"]:
        return
        # for g in geometry.geoms:
        #     show_geometry(g,ax=ax,show=False)
    else:
        print("no geometry/type to show")
        return
    if show:
        plt.show()


def save_animated_gif(image_files, output_filename, frame_duration_ms):
    # Zoek afbeeldingen in de huidige map

    if not image_files:
        print(
            "⚠️ Geen afbeeldingsbestanden gevonden in de map. Zorg dat er .png of .jpg bestanden aanwezig zijn."
        )
    else:
        frames = [Image.open(img) for img in image_files]
        # output_filename = "animated.gif"
        frame_duration_ms = 500  # duur per frame in milliseconden

        frames[0].save(
            output_filename,
            format="GIF",
            save_all=True,
            append_images=frames[1:],
            duration=frame_duration_ms,
            loop=0,
        )

    print(
        f"✅ Animated GIF succesvol opgeslagen als '{output_filename}' met {len(frames)} frames."
    )
