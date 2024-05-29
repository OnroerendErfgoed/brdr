import logging

import geopandas as gpd
import matplotlib.pyplot as plt


def show_results(results, results_diff_pos,results_diff_neg,thematic_dict, reference_dict):
    try:
        ax = gpd.GeoSeries(list(reference_dict.values())).plot(color='#FFF8C9',edgecolor='black',linewidth=2.0, zorder=1)
        ax_thematic_dict = gpd.GeoSeries(list(thematic_dict.values())).plot(ax=ax, alpha=0.8,color='none',hatch="/", edgecolor='#0000FF',linewidth=3.0,linestyle='dashdot', zorder=3)
        ax_result = gpd.GeoSeries(list(results.values())).plot(ax=ax,alpha=0.5,color='none',hatch=" ", edgecolor='green',linewidth=7.0, zorder=2)
        ax_diff_pos = gpd.GeoSeries(list(results_diff_pos.values())).plot(ax=ax, color='none', edgecolor='green', hatch="+",linewidth=0.0, linestyle='dashdot', zorder=4)
        ax_diff_neg = gpd.GeoSeries(list(results_diff_neg.values())).plot(ax=ax, color='none', edgecolor='red', hatch="+",linewidth=0.0, linestyle='dashdot', zorder=5)
        plt.show()
    except:
        logging.error("show_results: Error while showing results")
    return


def show_individual_results(results, results_diff_pos,results_diff_neg,thematic_dict, reference_dict,key):
    results = {key:results[key]}
    results_diff_pos = {key: results_diff_pos[key]}
    results_diff_neg = {key: results_diff_neg[key]}
    thematic_dict = {key: thematic_dict[key]}
    show_results(results, results_diff_pos,results_diff_neg,thematic_dict, reference_dict)
    return

def plot_diffs(series, results_diff):
    for key in results_diff:
        if len(results_diff[key]) == len(series):
            lst_diffs = list(results_diff[key].values())
            plt.plot(series, lst_diffs, label=str(key))
    plt.xlabel("relevant distance")
    plt.ylabel("difference")
    plt.title("Relevant distance vs difference")
    plt.legend()
    plt.show()
    return
