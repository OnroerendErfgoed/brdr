import os
import unittest

from shapely import Point
from shapely import from_wkt
from shapely.geometry import Polygon
from shapely.geometry import shape

from brdr.aligner import Aligner
from brdr.constants import AREA_ATTRIBUTE
from brdr.enums import GRBType, AlignerResultType
from brdr.enums import OpenDomainStrategy
from brdr.geometry_utils import (
    _grid_bounds,
    geom_from_wkt,
    geometric_equality,
    safe_equals,
)
from brdr.geometry_utils import buffer_neg_pos
from brdr.grb import (
    GRBActualLoader,
)
from brdr.loader import GeoJsonLoader, DictLoader
from brdr.typings import FeatureCollection, ProcessResult


class TestAligner(unittest.TestCase):

    def test_buffer_neg_pos(self):
        # Check if the result of the _buffer_neg_pos gives an equal geometry
        out = buffer_neg_pos(self.sample_geom, 3)
        self.assertEqual(self.sample_geom, out)

    def setUp(self):
        # Create a sample geometry for testing
        self.sample_aligner = Aligner()
        self.sample_geom = Polygon([(0, 0), (0, 10), (10, 10), (10, 0)])

    def test_buffer_neg_pos(self):
        # Check if the result of the _buffer_neg_pos gives an equal geometry
        out = buffer_neg_pos(self.sample_geom, 3)
        self.assertEqual(self.sample_geom, out)

    def test_grid_bounds_1(self):
        # Test _grid_bounds function
        delta = 1.0
        grid_partitions = _grid_bounds(self.sample_geom, delta)

        # Check if the result is a list of Polygon objects
        self.assertIsInstance(grid_partitions, list)
        for partition in grid_partitions:
            self.assertIsInstance(partition, Polygon)

        # Add more specific tests based on your requirements

    def test_grid_bounds_2(self):
        # Test _grid_bounds function
        delta = 2.0
        grid_partitions = _grid_bounds(self.sample_geom, delta)

        # Check if the result is a list of Polygon objects
        self.assertIsInstance(grid_partitions, list)
        for partition in grid_partitions:
            self.assertIsInstance(partition, Polygon)

    def test_export_results(self):
        aligner = Aligner()
        aligner.load_thematic_data(
            DictLoader(
                {"theme_id_1": from_wkt("POLYGON ((0 0, 0 9, 5 10, 10 0, 0 0))")}
            )
        )
        aligner.load_reference_data(
            DictLoader({"ref_id_1": from_wkt("POLYGON ((0 1, 0 10,8 10,10 1,0 1))")})
        )
        aligner.process()
        path = "./tmp/"
        resulttype = AlignerResultType.PROCESSRESULTS
        aligner.save_results(path=path, resulttype=resulttype)
        filenames = [
            resulttype.value + f"_{k}.geojson" for k in ProcessResult.__annotations__
        ]
        for file_name in os.listdir(path):
            os.remove(path + file_name)
            assert file_name in filenames
        os.rmdir(path)

    def test_get_formula_full_intersection(self):
        # Test when intersection equals reference geometry
        key = "a"
        ref_dict = {key: self.sample_geom}
        self.sample_aligner.load_reference_data(DictLoader(ref_dict))
        res = self.sample_aligner.get_brdr_formula(self.sample_geom, with_geom=True)
        self.assertTrue(res["full"])
        result = res["reference_features"][key]
        self.assertTrue(result["full"])
        self.assertEqual(result["percentage"], 100)

    def test_get_formula_partial_intersection(self):
        # Test when intersection is partial
        key = "a"
        ref_dict = {key: self.sample_geom.buffer(0.5)}
        self.sample_aligner.load_reference_data(DictLoader(ref_dict))
        res = self.sample_aligner.get_brdr_formula(self.sample_geom, with_geom=True)
        self.assertFalse(res["full"])
        result = res["reference_features"][key]
        self.assertFalse(result["full"])
        self.assertGreater(result["percentage"], 0)
        self.assertLess(result["percentage"], 100)

    def test_process_geometry(self):
        # Test if processed geometry is equal to reference geometry
        key_ref = "a"
        ref_dict = {key_ref: self.sample_geom}
        self.sample_aligner.load_reference_data(DictLoader(ref_dict))
        thematic_geom = self.sample_geom.buffer(0.5)
        process_result = self.sample_aligner.process_geometry(thematic_geom)
        self.assertTrue(
            from_wkt(process_result["result"].wkt).equals(
                from_wkt(self.sample_geom.wkt)
            )
        )
        assert geometric_equality(
            process_result.get("result"),
            self.sample_geom,
            self.sample_aligner.correction_distance,
            self.sample_aligner.mitre_limit,
        )

        self.assertFalse(process_result["result_diff"].is_empty)

    def test_line(self):
        aligner = Aligner(max_workers=-1)
        wkt = "MULTILINESTRING ((174024.1298775521281641 179420.42107488788315095, 174042.22504128722357564 179413.1830093938333448, 174044.36356063772109337 179418.28255553735652938, 174049.29860529277357273 179418.94056149138486944, 174050.28561422377242707 179422.72409572688047774, 174054.39815143629675731 179421.90158828438143246, 174057.52367971779312938 179420.75007786488276906, 174054.56265292479656637 179408.08346325031016022, 174047.16008594224695116 179398.54237691726302728, 174036.960993655200582 179404.29992901478544809, 174032.84845644267625175 179400.51639477928983979, 174029.88742964965058491 179397.88437096326379105))"
        thematic_dict = {"theme_id_1": from_wkt(wkt)}

        loader = DictLoader(thematic_dict)
        aligner.load_thematic_data(loader)
        loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
        aligner.load_reference_data(loader)
        relevant_distance = 3
        dict_processresults = aligner.process(relevant_distance=relevant_distance)
        self.assertEqual(
            dict_processresults["theme_id_1"][relevant_distance]["result"].geom_type,
            "LineString",
        )

    def test_reference_mix(self):
        "reference exists out of points, lines and polygons"
        aligner = Aligner(max_workers=-1)
        wkt = "POLYGON ((174043.1035556931165047 179299.26716663906699978, 174042.53709723605425097 179306.22651339753065258, 174048.52537235378986225 179305.9837454873486422, 174055.16102856534416787 179305.49820966698462144, 174055.16102856534416787 179293.19796888460405171, 174071.26463327385135926 179293.92627261512097903, 174071.34555591057869606 179287.69522958720335737, 174083.40302878280635923 179282.83987138362135738, 174072.154782277852064 179278.14635845352313481, 174072.72124073494342156 179260.58614628392388113, 174056.69855866313446313 179258.88677091267891228, 174048.60629499051719904 179273.21007761321379803, 174038.97650122008053586 179273.85745870703249238, 174041.97063877896289341 179286.72415794650441967, 174044.96477633781614713 179292.71243306424003094, 174035.25405993068125099 179292.14597460714867339, 174043.1035556931165047 179299.26716663906699978))"
        thematic_dict = {"theme_id_1": from_wkt(wkt)}
        loader = DictLoader(thematic_dict)
        aligner.load_thematic_data(loader)
        point_1 = geom_from_wkt(
            "POINT (174043.75093678693519905 179298.98393741055042483)"
        )
        point_2 = geom_from_wkt(
            "POINT (174079.35689694649772719 179283.04217797549790703)"
        )
        line_1 = geom_from_wkt(
            "LINESTRING (174021.41628905048128217 179291.53905483175185509, 174081.29904022792470641 179294.45226975387777202)"
        )
        line_2 = geom_from_wkt(
            "LINESTRING (174039.13834649353520945 179282.79941006531589665, 174038.16727485283627175 179272.92684838469722308, 174052.65242682682583109 179273.33146156833390705, 174093.43743573685060255 179285.46985707725980319)"
        )
        polygon_1 = geom_from_wkt(
            "POLYGON ((174029.67039799655321985 179306.67158789956010878, 174048.44444971703342162 179306.18605207919608802, 174048.92998553739744239 179258.92723223107168451, 174030.15593381691724062 179259.41276805143570527, 174029.67039799655321985 179306.67158789956010878))"
        )
        polygon_2 = geom_from_wkt(
            "POLYGON ((174054.91826065513305366 179305.53867098540649749, 174056.37486811622511595 179258.76538695761701092, 174072.72124073491431773 179260.38383969213464297, 174070.61725218003266491 179306.99527844646945596, 174054.91826065513305366 179305.53867098540649749))"
        )
        reference_dict = {
            "ref_id_1": point_1,
            "ref_id_2": point_2,
            "ref_id_3": line_1,
            "ref_id_4": line_2,
            "ref_id_5": polygon_1,
            "ref_id_6": polygon_2,
        }

        loader = DictLoader(reference_dict)
        aligner.load_reference_data(loader)
        relevant_distance = 5
        dict_processresults = aligner.process(relevant_distance=relevant_distance)
        # TODO-check!
        # self.assertEqual(
        #     dict_processresults["theme_id_1"][relevant_distance]["result"].geom_type,
        #     "GeometryCollection",
        # )

    def test_load_reference_data_grb_actual_adp(self):
        thematic_dict = {
            "theme_id_1": from_wkt(
                "MultiPolygon (((174184.09476602054201066 171899.68933439542888664, "
                "174400.56834639035514556 171832.959863749332726, "
                "174388.65236948925303295 171770.99678386366576888, "
                "174182.10876987033407204 171836.13745758961886168, "
                "174184.88916448061354458 171873.07698598300339654, "
                "174184.09476602054201066 171899.68933439542888664)))"
            )
        }
        self.sample_aligner.load_thematic_data(DictLoader(thematic_dict))
        # LOAD REFERENCE DICTIONARY
        self.sample_aligner.load_reference_data(
            GRBActualLoader(
                aligner=self.sample_aligner, grb_type=GRBType.ADP, partition=1000
            )
        )
        self.assertGreater(len(self.sample_aligner.dict_reference), 0)

    def test_load_reference_data_grb_actual_gbg(self):
        thematic_dict = {
            "theme_id_1": from_wkt(
                "MultiPolygon (((174184.09476602054201066 171899.68933439542888664, "
                "174400.56834639035514556 171832.959863749332726, "
                "174388.65236948925303295 171770.99678386366576888, "
                "174182.10876987033407204 171836.13745758961886168, "
                "174184.88916448061354458 171873.07698598300339654, "
                "174184.09476602054201066 171899.68933439542888664)))"
            )
        }
        self.sample_aligner.load_thematic_data(DictLoader(thematic_dict))
        # LOAD REFERENCE DICTIONARY
        self.sample_aligner.load_reference_data(
            GRBActualLoader(
                aligner=self.sample_aligner, grb_type=GRBType.GBG, partition=1000
            )
        )
        self.assertGreater(len(self.sample_aligner.dict_reference), 0)

    def test_load_reference_data_grb_actual_knw(self):
        thematic_dict = {
            "theme_id_1": from_wkt(
                "MultiPolygon (((174184.09476602054201066 171899.68933439542888664, "
                "174400.56834639035514556 171832.959863749332726, "
                "174388.65236948925303295 171770.99678386366576888, "
                "174182.10876987033407204 171836.13745758961886168, "
                "174184.88916448061354458 171873.07698598300339654, "
                "174184.09476602054201066 171899.68933439542888664)))"
            )
        }
        self.sample_aligner.load_thematic_data(DictLoader(thematic_dict))
        # LOAD REFERENCE DICTIONARY
        self.sample_aligner.load_reference_data(
            GRBActualLoader(
                aligner=self.sample_aligner, grb_type=GRBType.KNW, partition=1000
            )
        )
        self.sample_aligner.process()
        self.assertGreaterEqual(len(self.sample_aligner.dict_reference), 0)

    def test_all_od_strategies(self):
        # Load thematic data & reference data
        thematic_dict = {
            "theme_id_1": from_wkt("POLYGON ((0 0, 0 9, 5 10, 10 0, 0 0))")
        }
        # ADD A REFERENCE POLYGON TO REFERENCE DICTIONARY
        reference_dict = {"ref_id_1": from_wkt("POLYGON ((0 1, 0 10,8 10,10 1,0 1))")}
        # LOAD THEMATIC DICTIONARY
        self.sample_aligner.load_thematic_data(DictLoader(thematic_dict))
        # LOAD REFERENCE DICTIONARY
        self.sample_aligner.load_reference_data(DictLoader(reference_dict))
        relevant_distance = 1
        for od_strategy in OpenDomainStrategy:
            process_result = self.sample_aligner.process(
                relevant_distance=relevant_distance,
                od_strategy=od_strategy,
                threshold_overlap_percentage=50,
            )
            self.assertEqual(len(process_result["theme_id_1"][relevant_distance]), 7)

    def test_process_interior_ring(self):
        thematic_dict = {
            "theme_id_1": from_wkt(
                "MultiPolygon (((174370.67910432978533208 171012.2469546866195742, "
                "174461.24052877808571793 171002.71417316573206335, "
                "174429.46459037516615354 170869.25523187351063825, "
                "174373.85669817010057159 170894.6759825958579313, "
                "174369.09030740967136808 170896.6619787460367661, "
                "174370.67910432978533208 171012.2469546866195742),"
                "(174400.07184735251939856 170963.78864862219779752, "
                "174401.26344504262669943 170926.4519209987774957, "
                "174419.53460962430108339 170922.47992869839072227, "
                "174429.06739114515949041 170958.62505863170372322, "
                "174400.07184735251939856 170963.78864862219779752)))"
            )
        }
        self.sample_aligner.load_thematic_data(DictLoader(thematic_dict))
        # LOAD REFERENCE DICTIONARY
        self.sample_aligner.load_reference_data(
            GRBActualLoader(
                aligner=self.sample_aligner, grb_type=GRBType.GBG, partition=1000
            )
        )
        result_dict = self.sample_aligner.process()
        self.assertEqual(len(result_dict), len(thematic_dict))

    def test_process_no_added_area(self):
        aligner = Aligner(crs="EPSG:31370", multi_as_single_modus=True)
        # ADD A THEMATIC POLYGON TO THEMATIC DICTIONARY and LOAD into Aligner
        id = "my_theme_id"
        wkt = "MULTIPOLYGON (((141894.30860000103712082 192030.99359999969601631, 141895.30039999634027481 192030.43499999865889549, 141916.94680000096559525 192018.24240000173449516, 141935.11959999799728394 192008.00629999861121178, 141938.46660000085830688 192006.36679999902844429, 141938.01630000025033951 192004.09310000017285347, 141936.13579999655485153 191997.07739999890327454, 141927.83190000057220459 191974.37330000102519989, 141924.48290000110864639 191964.26909999921917915, 141918.48789999634027481 191946.18219999969005585, 141897.54580000042915344 191892.38010000064969063, 141885.81229999661445618 191863.05810000002384186, 141883.73430000245571136 191858.32620000094175339, 141881.57930000126361847 191851.27149999886751175, 141880.02970000356435776 191844.1534000001847744, 141875.37430000305175781 191822.76999999955296516, 141875.30219999700784683 191822.43849999830126762, 141875.31379999965429306 191822.42839999869465828, 141875.32599999755620956 191822.41780000180006027, 141873.66679999977350235 191810.330400001257658, 141873.53939999639987946 191810.47140000015497208, 141824.16929999738931656 191865.09059999883174896, 141807.49130000174045563 191883.54179999977350235, 141806.4221000000834465 191884.39990000054240227, 141804.12879999727010727 191884.68789999932050705, 141802.18280000239610672 191884.19500000029802322, 141799.86169999837875366 191882.65430000051856041, 141797.43410000205039978 191879.70760000124573708, 141796.65850000083446503 191879.23550000041723251, 141795.6005999967455864 191879.14139999821782112, 141793.61129999905824661 191879.65769999846816063, 141791.29720000177621841 191881.41699999943375587, 141775.49679999798536301 191898.58630000054836273, 141762.63560000061988831 191912.02459999918937683, 141746.34480000287294388 191931.99890000000596046, 141739.71980000287294388 191941.33379999920725822, 141718.17260000109672546 191963.87099999934434891, 141707.33889999985694885 191974.1517999991774559, 141691.99570000171661377 191990.56419999897480011, 141682.7752000018954277 191999.544599998742342, 141680.59319999814033508 192001.8539000004529953, 141677.57880000025033951 192009.98910000175237656, 141677.57880000025033951 192024.72419999912381172, 141687.39810000360012054 192032.38650000095367432, 141710.51529999822378159 192048.2443000003695488, 141718.70740000158548355 192054.53130000084638596, 141792.93550000339746475 192111.49830000102519989, 141790.39050000160932541 192115.83139999955892563, 141771.93940000236034393 192145.56280000135302544, 141756.72280000150203705 192171.94390000030398369, 141750.73979999870061874 192185.24599999934434891, 141746.87189999967813492 192192.59640000015497208, 141744.20000000298023224 192197.03090000152587891, 141741.65730000287294388 192203.05490000173449516, 141738.98480000346899033 192211.15890000015497208, 141736.27549999952316284 192222.19130000099539757, 141734.4122999981045723 192232.36960000172257423, 141734.01650000363588333 192233.20839999988675117, 141733.4156000018119812 192233.86070000007748604, 141732.60769999772310257 192234.36529999971389771, 141731.69009999930858612 192234.72720000147819519, 141730.37650000303983688 192234.86820000037550926, 141727.74570000171661377 192234.72109999880194664, 141708.45470000058412552 192234.2791999988257885, 141680.76659999787807465 192234.02210000157356262, 141656.07289999723434448 192234.38520000129938126, 141629.07209999859333038 192233.29850000143051147, 141618.51850000023841858 192231.88439999893307686, 141603.34749999642372131 192228.10990000143647194, 141597.26340000331401825 192226.65670000016689301, 141592.5292000025510788 192225.30229999870061874, 141568.39479999989271164 192219.11100000143051147, 141554.37330000102519989 192214.54509999975562096, 141541.26389999687671661 192208.97210000082850456, 141529.27139999717473984 192202.5982000008225441, 141527.42339999973773956 192207.54980000108480453, 141526.95019999891519547 192208.81769999861717224, 141538.93559999763965607 192213.12829999998211861, 141551.88490000367164612 192218.77100000157952309, 141563.40359999984502792 192223.01000000163912773, 141568.94730000197887421 192224.56370000168681145, 141588.22190000116825104 192230.68310000002384186, 141599.65079999715089798 192233.4813000001013279, 141604.27549999952316284 192234.52360000088810921, 141612.59619999676942825 192236.39869999885559082, 141625.6671999990940094 192238.60940000042319298, 141640.7581000030040741 192239.92980000004172325, 141687.70180000364780426 192239.60280000045895576, 141708.04810000211000443 192239.78669999912381172, 141725.32699999958276749 192239.62020000070333481, 141736.44510000199079514 192239.40599999949336052, 141738.10999999940395355 192233.08529999852180481, 141744.46729999780654907 192208.94959999993443489, 141747.24490000307559967 192201.90359999984502792, 141750.76910000294446945 192194.79340000078082085, 141757.62960000336170197 192183.66059999912977219, 141765.30849999934434891 192168.3214000016450882, 141766.5877000018954277 192166.07369999960064888, 141779.35530000180006027 192173.02180000022053719, 141801.67239999771118164 192199.32169999927282333, 141831.40259999781847 192121.48229999840259552, 141831.68190000206232071 192121.63430000096559525, 141839.40439999848604202 192125.83700000122189522, 141847.05669999867677689 192105.80189999938011169, 141848.41929999738931656 192102.23440000042319298, 141852.30349999666213989 192092.06469999998807907, 141856.9464000016450882 192074.37240000069141388, 141855.61249999701976776 192070.44669999927282333, 141854.58449999988079071 192067.42040000110864639, 141851.04050000011920929 192056.98910000175237656, 141853.14169999957084656 192055.68459999933838844, 141855.8546999990940094 192053.98660000041127205, 141857.53270000219345093 192052.94280000030994415, 141865.36980000138282776 192048.06769999861717224, 141878.53230000287294388 192039.87979999929666519, 141877.30830000340938568 192037.69469999894499779, 141888.70970000326633453 192031.27279999852180481, 141889.19900000095367432 192031.21249999850988388, 141890.32779999822378159 192033.23580000177025795, 141894.30860000103712082 192030.99359999969601631),(141798.70589999854564667 192111.06799999997019768, 141801.00859999656677246 192106.89739999920129776, 141829.50010000169277191 192122.21880000084638596, 141800.99939999729394913 192196.83920000120997429, 141781.65139999985694885 192173.08870000019669533, 141767.0609000027179718 192165.24260000139474869, 141773.47720000147819519 192153.96849999949336052, 141785.42890000343322754 192135.11470000073313713, 141798.70589999854564667 192111.06799999997019768)))"
        thematic_dict = {id: geom_from_wkt(wkt)}
        loader = DictLoader(data_dict=thematic_dict)
        aligner.load_thematic_data(loader)

        # Load reference data: The actual GRB-parcels
        loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
        aligner.load_reference_data(loader)
        # EXECUTE THE ALIGNMENT
        relevant_distance = 5
        process_result = aligner.process(
            relevant_distance=relevant_distance,
            od_strategy=OpenDomainStrategy.SNAP_INNER_SIDE,
            threshold_overlap_percentage=50,
        )

        result = process_result[id][relevant_distance]["result"]
        result_diff_plus = process_result[id][relevant_distance]["result_diff_plus"]
        result_diff_min = process_result[id][relevant_distance]["result_diff_min"]
        assert not result.is_empty
        assert result_diff_plus.is_empty
        assert not result_diff_min.is_empty

    def test_process_big_relevant_distance(self):
        wkt = "POLYGON ((172791.07652649749070406 171319.35066181665752083, 172821.75881976058008149 171308.36308382378774695, 172838.13653035368770361 171321.42378973981249146, 172879.59908881731098518 171303.80220239277696237, 172842.49009899239172228 171205.95056441865745001, 172804.55185799815808423 171219.42589591932483017, 172789.21071136664249934 171221.29171105020213872, 172757.69916693426785059 171270.63215562188997865, 172766.4063042116467841 171279.13198010693304241, 172776.56463103523128666 171296.96088024630444124, 172785.63235233756131493 171305.36084497725823894, 172782.68072999455034733 171312.79389392648590729, 172791.07652649749070406 171319.35066181665752083))"
        thematic_dict = {"theme_id_1": from_wkt(wkt)}
        self.sample_aligner.load_thematic_data(DictLoader(thematic_dict))
        # LOAD REFERENCE DICTIONARY
        self.sample_aligner.load_reference_data(
            GRBActualLoader(
                aligner=self.sample_aligner, grb_type=GRBType.ADP, partition=1000
            )
        )
        result_dict = self.sample_aligner.process(
            od_strategy=OpenDomainStrategy.SNAP_ALL_SIDE, relevant_distance=100
        )
        self.assertEqual(len(result_dict), len(thematic_dict))

    def test_process_circle(self):
        geometry = Point(0, 0).buffer(3)
        # geometry = MultiPolygon([geometry])
        thematic_dict = {"key": geometry}
        self.sample_aligner.load_thematic_data(DictLoader(thematic_dict))
        # LOAD REFERENCE DICTIONARY
        self.sample_aligner.load_reference_data(
            GRBActualLoader(
                aligner=self.sample_aligner, grb_type=GRBType.GBG, partition=1000
            )
        )
        relevant_distance = 1
        results_dict = self.sample_aligner.process(relevant_distance=relevant_distance)
        self.assertEqual(geometry, results_dict["key"][relevant_distance]["result"])

    def test__prepare_thematic_data(self):
        aligner = Aligner()
        geojson: FeatureCollection = {
            "type": "FeatureCollection",
            "name": "theme",
            "crs": {
                "type": "name",
                "properties": {"name": "urn:ogc:def:crs:EPSG::31370"},
            },
            "features": [
                {
                    "type": "Feature",
                    "properties": {"fid": 4, "id": 4, "theme_identifier": "4"},
                    "geometry": {
                        "type": "MultiPolygon",
                        "coordinates": [
                            [
                                [
                                    [174647.924166895216331, 170964.980246312363306],
                                    [174693.204879119351972, 170943.531487890402786],
                                    [174696.382472959638108, 170930.82111252922914],
                                    [174678.905706838035258, 170901.428369506524177],
                                    [174660.634542256389977, 170861.708446502889274],
                                    [174641.56897921464406, 170820.399726579111302],
                                    [174593.905071610264713, 170839.46528962085722],
                                    [174614.559431572153699, 170881.568408004706725],
                                    [174628.064205393398879, 170926.054721768799936],
                                    [174647.924166895216331, 170964.980246312363306],
                                ]
                            ]
                        ],
                    },
                }
            ],
        }
        thematic_loader = GeoJsonLoader(_input=geojson, id_property="theme_identifier")
        aligner.dict_thematic, properties, source = thematic_loader.load_data()
        assert aligner.dict_thematic == {"4": shape(geojson["features"][0]["geometry"])}
        self.assertGreater(len(aligner.dict_thematic), 0)

    def test_get_reference_as_geojson(self):
        self.sample_aligner.load_thematic_data(
            DictLoader(
                {"theme_id_1": from_wkt("POLYGON ((0 0, 0 9, 5 10, 10 0, 0 0))")}
            )
        )
        self.sample_aligner.load_reference_data(
            DictLoader({"ref_id_1": from_wkt("POLYGON ((0 1, 0 10,8 10,10 1,0 1))")})
        )
        self.sample_aligner.process()
        self.sample_aligner.get_input_as_geojson()

    def test_fully_aligned_input(self):
        aligned_shape = from_wkt("POLYGON ((0 0, 0 9, 5 10, 10 0, 0 0))")
        self.sample_aligner.load_thematic_data(
            DictLoader({"theme_id_1": aligned_shape})
        )
        self.sample_aligner.load_reference_data(DictLoader({"ref_id_1": aligned_shape}))
        relevant_distance = 1
        result = self.sample_aligner.process(relevant_distance=relevant_distance)
        assert safe_equals(
            result["theme_id_1"][relevant_distance].get("result"), aligned_shape
        )
        assert geometric_equality(
            result["theme_id_1"][relevant_distance].get("result"),
            aligned_shape,
            self.sample_aligner.correction_distance,
            self.sample_aligner.mitre_limit,
        )
        assert result["theme_id_1"][relevant_distance].get("result_diff").is_empty
        assert result["theme_id_1"][relevant_distance].get("result_diff_min").is_empty
        assert result["theme_id_1"][relevant_distance].get("result_diff_plus").is_empty

    def test_remark_for_poly_multipoly(self):
        shape = from_wkt(
            "MultiPolygon(((48893.03662109375 214362.93756103515625, 48890.8258056640625 214368.482666015625, 48890.7159423828125 214368.44110107421875, 48887.6488037109375 214367.2845458984375, 48886.3800048828125 214368.68017578125, 48885.1068115234375 214370.08062744140625, 48884.3330078125 214369.782470703125, 48882.563720703125 214369.10064697265625, 48882.1116943359375 214370.1346435546875, 48878.5626220703125 214368.70196533203125, 48877.839111328125 214368.40997314453125, 48877.2352294921875 214369.79376220703125, 48876.7911376953125 214369.60687255859375, 48875.0850830078125 214373.62353515625, 48875.478759765625 214373.8182373046875, 48881.5286865234375 214376.81109619140625, 48885.10546875 214372.36151123046875, 48887.0050048828125 214370.08538818359375, 48888.4698486328125 214368.330078125, 48890.366943359375 214369.2685546875, 48901.0638427734375 214374.56024169921875, 48905.0159912109375 214369.61175537109375, 48904.472900390625 214367.53851318359375, 48893.03662109375 214362.93756103515625)))"
        )
        self.sample_aligner.load_thematic_data(DictLoader({"theme_id_1": shape}))
        self.sample_aligner.load_reference_data(
            GRBActualLoader(
                grb_type=GRBType.ADP, partition=1000, aligner=self.sample_aligner
            )
        )
        rd = 2
        self.sample_aligner.process(relevant_distances=[rd])
        assert (
            self.sample_aligner.dict_processresults["theme_id_1"][rd]["remark"]
            == " | Difference in amount of geometries"
        )

    def test_fully_aligned_geojson_output(self):
        aligned_shape = from_wkt(
            "MultiPolygon (((173463.11530961000244133 174423.83310307000647299, "
            "173460.22633100001257844 174422.02316300000529736, "
            "173455.24681099998997524 174429.98009100000490434, "
            "173454.4299790000077337 174429.34482699999352917, "
            "173452.06690700000035577 174432.43058700000983663, "
            "173451.25743500000680797 174431.8672589999914635, "
            "173448.74844299998949282 174434.96249100001296028, "
            "173448.5809550000121817 174435.80485899999621324, "
            "173455.39772300000186078 174441.47852299999794923, "
            "173461.44169100001454353 174446.50898700000834651, "
            "173472.15932299999985844 174429.49919500001124106, "
            "173466.18524341000011191 174425.75641125999391079, "
            "173466.18524347001221031 174425.75641117000486702, "
            "173463.11530969000887126 174423.83310300001176074, "
            "173463.11530961000244133 174423.83310307000647299)))"
        )

        self.sample_aligner.load_thematic_data(
            DictLoader({"theme_id_1": aligned_shape})
        )
        self.sample_aligner.load_reference_data(DictLoader({"ref_id_1": aligned_shape}))
        self.sample_aligner.process()
        fcs = self.sample_aligner.get_results_as_geojson(formula=True)
        assert fcs["result"]["features"][0]["properties"][AREA_ATTRIBUTE] > 0
        assert fcs["result_diff"]["features"][0]["properties"][AREA_ATTRIBUTE] == 0
        assert fcs["result_diff_min"]["features"][0]["properties"][AREA_ATTRIBUTE] == 0
        assert fcs["result_diff_plus"]["features"][0]["properties"][AREA_ATTRIBUTE] == 0
