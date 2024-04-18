import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("../tests/output/stats1710857407.9886606.csv", sep=";")
pd.options.display.max_rows = 5
print(pd.options.display.max_rows)
print(df)

df["big_area"] = df["area"] > 10000  # = 10 000m² is 2 voetbalvelden
# df = df.assign(max = df.groupby(['key','strategy','full_percentage'])[
# 'diff'].transform('max')).values
df["max"] = df.groupby(["key", "strategy", "full_percentage"])["diff"].transform("max")
print(df)
test = df.describe()
print(test)
# profile = ProfileReport(df)
# profile.to_file('test.html')

result = df
# result = result.loc[(df['distance'] != 10) & (df['strategy'] == 1) & (df[
# 'full_percentage'] == 50)]
result = result.loc[(df["area"] < 10000)]
result = result.loc[df["diff"] == df["max"]]
print(result)
result = result.sort_values(
    ["distance", "key", "strategy", "full_percentage", "big_area", "max"]
).drop_duplicates(["key", "strategy", "full_percentage", "big_area", "max"])
# result =result.loc[result.groupby(['key','strategy','full_percentage','max'])[
# 'distance'].idxmin()]

print(result)

plt.scatter(result["max"], result["area"] / 1000)
plt.show()
# profile = ProfileReport(result)
# profile.to_file('test.html')

ax = result.hist(
    column="distance",
    by=["strategy", "full_percentage", "big_area"],
    bins=50,
    grid=False,
    figsize=(16, 20),
    layout=(2, 1),
    sharex=True,
    color="#86bf91",
    zorder=2,
    rwidth=0.9,
)

# for i,x in enumerate(ax):
#     # Despine
#     # x.spines['right'].set_visible(False)
#     # x.spines['top'].set_visible(False)
#     # x.spines['left'].set_visible(False)
#
#     # Switch off ticks

#     #x.tick_params(
#       axis="both",
#       which="both",
#       bottom="off",
#       top="off",
#       labelbottom="on",
#       left="off",
#       right="off",
#       labelleft="on"
#       )
#
#     # # Draw horizontal axis lines
#     # vals = x.get_yticks()
#     # for tick in vals:
#     #     x.axhline(y=tick, linestyle='dashed', alpha=0.4, color='#eeeeee', zorder=1)
#
#     # Set x-axis label
#     x.set_xlabel("Distance", labelpad=20, weight='bold', size=12)
#
#     # Set y-axis label
#     if i == 1:
#         x.set_ylabel("Count", labelpad=50, weight='bold', size=12)
#
#     #  Format y-axis label
#     #x.yaxis.set_major_formatter(StrMethodFormatter('{x:,g}'))
#
#     #x.tick_params(axis='x', rotation=0)
# #
plt.show()

result_bplt = result.boxplot(
    column="distance", by=["strategy", "full_percentage", "big_area"]
)

result_bplt.plot()
plt.show()

# CODE to create stats.csv
# time = str(time.time())
# with open('stats' + time + '.csv', 'w', newline='') as csvfile:
#     writer = csv.writer(csvfile, delimiter=';')
# ,quotechar='"', quoting=csv.QUOTE_MINIMAL)
#     writer.writerow(
#     ['distance',
#     'strategy',
#     'full_percentage',
#     'key','area',
#     'diff_plus',
#     'diff_min',
#     'diff','diff_prc']
#     )
#     for full_percentage in [0,25,50,75,100]:
#         logging.info('full percentage: ' + str(full_percentage))
#         for od in [-1, 0, 1, 2]:
#             logging.info('od_strategy: ' + str(od))
#             for s in [0.2,0.5,1,1.5,2,3,4,5,6,8,10]:
#                 logging.info('relevant_distance: ' + str(s))
#                 results,results_diff,results_diff_plus,results_diff_min = x.process(
#                 s,od,full_percentage
#                 )
#                 for key in results:
#                     results_area = results[key].area
#                     results_diff_plus_area = results_diff_plus[key].area
#                     results_diff_min_area = results_diff_min[key].area
#                     diff_area = results_diff_plus_area + results_diff_min_area
#                     diff_area_prc = diff_area*100/results_area
#                      writer.writerow(
#                             [
#                                 s,
#                                 od,
#                                 full_percentage,
#                                 key,
#                                 results_area,
#                                 results_diff_plus_area,
#                                 results_diff_min_area,
#                                 diff_area,
#                                 diff_area_prc,
#                             ]
#                         )
