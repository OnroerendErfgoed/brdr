import geopandas as gpd
import matplotlib.pyplot as plt


def show_results(results, results_diff_pos,results_diff_neg):
    results_array = list(results.values())
    results_diff_pos_array = list(results_diff_pos.values())
    results_diff_neg_array = list(results_diff_neg.values())

    results_geoseries = gpd.GeoSeries(results_array)
    results_diff_pos_geoseries = gpd.GeoSeries(results_diff_pos_array)
    ax1 = results_diff_pos_geoseries.plot(color='none', edgecolor='green', hatch="/", zorder=3)
    results_diff_neg_geoseries = gpd.GeoSeries(results_diff_neg_array)
    ax2 = results_diff_neg_geoseries.plot(ax =ax1,color='none', edgecolor='red', hatch="/", zorder=2)
    results_geoseries.plot(ax=ax2, zorder=1)
    plt.show()
    return


def plot_diffs(series, results_diff):
    for key in results_diff:
        if len(results_diff[key]) == len(series):
            lst_diffs = list(results_diff[key].values())
            plt.plot(series, lst_diffs, label=str(key))
    plt.xlabel("relevante afstand")
    plt.ylabel("diff")
    plt.title("relevante afstand vs diff")
    plt.legend()
    plt.show()
    return
