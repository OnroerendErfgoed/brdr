import statistics
from datetime import datetime

from brdr.aligner import Aligner
from brdr.enums import OpenDomainStrategy, GRBType
from brdr.geometry_utils import geom_from_wkt
from brdr.grb import GRBActualLoader
from brdr.loader import GeoJsonFileLoader, DictLoader


def main():
    """
    EXAMPLE of a test to measure the speed of the aligner
    :return:
    """
    # Initiate brdr
    aligner = Aligner(max_workers=None)
    iterations = 2
    od_strategy = OpenDomainStrategy.SNAP_PREFER_VERTICES
    relevant_distance = 3
    grb_loader = True
    aligner.multi_as_single_modus = True
    # Load local thematic data and reference data
    # loader = GeoJsonFileLoader(
    #     "../tests/testdata/theme.geojson", "theme_identifier"
    # )
    # loader = GeoJsonFileLoader(
    #     "../tests/testdata/themelayer_not_referenced.geojson", "theme_identifier"
    # )
    wkt = "Polygon ((174159.86498764005955309 179307.03392354605603032, 174261.82966211286839098 179298.04720647388603538, 174257.68194654109538533 179230.99247139683575369, 174246.62137168302433565 179223.0426832175871823, 174246.27572871872689575 179216.47546689561568201, 174258.7188754340459127 179214.40160910971462727, 174248.69522946892539039 179177.07216896372847259, 174229.68486643160576932 179182.60245639277854934, 174198.92264260759111494 179199.88460460852365941, 174153.9890572466829326 179211.63646539521869272, 174141.20026756706647575 179215.09289503836771473, 174141.20026756706647575 179215.09289503836771473, 174148.45876981766195968 179252.76797814865130931, 174159.86498764005955309 179307.03392354605603032))"
    # big area
    # grb_loader = False
    # wkt = "Polygon((173773.66988192888675258 179514.52399982791393995, 173801.10608616709941998 179530.65382982976734638, 173931.40736798732541502 179611.27281510829925537, 174078.04176878236467019 179629.09588310681283474, 174280.05469680932583287 179569.44558737333863974, 174351.46164843259612098 179496.61523636244237423, 174335.60671524403733201 179455.38617126829922199, 174331.16058107223943807 179425.64216128829866648, 174348.64269148261519149 179363.99040084332227707, 174360.42594051550258882 179295.45192713756114244, 174341.39866621946566738 179197.04401163011789322, 174358.8869952189270407 179134.24635295104235411, 174283.39069474133430049 179099.59578586928546429, 174289.62856128191924654 178991.01253041252493858, 174167.19324753561522812 178942.42442009877413511, 174006.65193600184284151 178953.11535596195608377, 173714.49524312707944773 178916.33480669092386961, 173711.43248903870698996 179085.52870534732937813, 173728.93523737252689898 179269.6728869853541255, 173674.80573671424644999 179325.4449802627786994, 173607.21246577176498249 179317.13567749317735434, 173596.70450016154791228 179359.38730681408196688, 173674.31222398031968623 179434.04803174268454313, 173773.66988192888675258 179514.52399982791393995))"
    # wkt = "Polygon ((172594.16427666315576062 177732.38409853109624237, 176756.43182838105713017 177257.73955315977218561, 176464.34287738331477158 172310.48294563542003743, 171316.27511604799656197 172675.59413438258343376, 171316.27511604799656197 172675.59413438258343376, 171772.66410198199446313 174446.38339980645105243, 172594.16427666315576062 177732.38409853109624237))"
    geom = geom_from_wkt(wkt)
    thematic_dict = {"theme_id_1": geom}
    loader = DictLoader(thematic_dict)
    aligner.load_thematic_data(loader)
    if not grb_loader:
        loader = GeoJsonFileLoader(
            "../tests/testdata/reference_leuven.geojson", "capakey"
        )
    else:
        loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
    aligner.load_reference_data(loader)
    aligner.partial_snapping = True
    aligner.process(relevant_distance=relevant_distance, od_strategy=od_strategy)
    times = []
    total_starttime = datetime.now()
    for iter in range(1, iterations + 1):
        starttime = datetime.now()

        # Example how to use the Aligner
        aligner.process(relevant_distance=relevant_distance, od_strategy=od_strategy)
        # fcs = aligner.get_results_as_geojson(formula=True)
        endtime = datetime.now()
        seconds = (endtime - starttime).total_seconds()
        times.append(seconds)
        print(seconds)
    total_endtime = datetime.now()
    total_seconds = (total_endtime - total_starttime).total_seconds()
    print("Total time: " + str(total_seconds))
    print("duration: " + str(times))
    print("Min: " + str(min(times)))
    print("Max: " + str(max(times)))
    print("Mean: " + str(statistics.mean(times)))
    print("Median: " + str(statistics.median(times)))
    print("Stdv: " + str(statistics.stdev(times)))


# #BEFORE REFACTORING dict_series
# duration: [25.652311, 27.894154, 19.641618, 19.929254, 44.754033, 25.218422, 23.167992, 18.649832, 22.899336, 52.108296]
# Min: 18.649832
# Max: 52.108296
# Mean: 27.9915248
# Median: 24.193207
# Stdv: 11.28891821173264

# #AFTER refactoring
# duration: [21.313991, 16.558168, 16.590126, 18.111118, 16.872433, 17.928071, 18.32295, 17.87116, 19.516652, 16.729241]
# Min: 16.558168
# Max: 21.313991
# Mean: 17.981391
# Median: 17.8996155
# Stdv: 1.504459449440969

if __name__ == "__main__":
    main()
