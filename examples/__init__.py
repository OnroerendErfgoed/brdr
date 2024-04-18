import geopandas as gpd
import matplotlib.pyplot as plt


def show_results(results, results_diff):
    results_array = list(results.values())
    results_diff_array = list(results_diff.values())

    results_geoseries = gpd.GeoSeries(results_array)
    results_diff_geoseries = gpd.GeoSeries(results_diff_array)
    ax = results_diff_geoseries.plot(color="k", zorder=2)
    results_geoseries.plot(ax=ax, zorder=1)
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
