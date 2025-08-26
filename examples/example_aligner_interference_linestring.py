from brdr.aligner import Aligner
from brdr.enums import GRBType
from brdr.geometry_utils import geom_from_wkt
from brdr.grb import GRBActualLoader
from brdr.loader import DictLoader

# Press the green button in the gutter to run the script.
if __name__ == "__main__":
    """
    EXAMPLE to check interference when Aligner-class is re-used
    """
    # Initiate an Aligner
    aligner = Aligner(max_workers=-1)
    # Load thematic data & reference data

    wkt = "LINESTRING (171741.11190000033820979 171839.01070547936251387, 171751.68948904142598622 171847.23771917796693742, 171762.26707808251376264 171855.69979041084297933, 171772.72713835645117797 171862.9865739724773448, 171783.53978493178146891 171870.8610013697471004, 171794.94007534274714999 171877.79519862998859026, 171801.87427260301774368 171884.2592808217741549)"
    thematic_dict = {"id_1": geom_from_wkt(wkt)}
    loader = DictLoader(data_dict=thematic_dict)

    aligner.load_thematic_data(loader)
    # Load reference data
    loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
    aligner.load_reference_data(loader)

    # PREDICT the 'stable' relevant distances, for a series of relevant distances

    # predict which relevant distances are interesting to propose as resulting geometry
    dict_evaluated1, prop_dictionary1 = aligner.evaluate()

    wkt = "MULTILINESTRING ((174135.22254434687783942 179531.39341345359571278, 174121.64978590866667219 179537.66386910632718354, 174118.0780073722708039 179538.69571623904630542, 174121.25292162684490904 179549.09356042271247134, 174114.98246597408433445 179551.23662754453835078, 174116.56992310137138702 179556.95147320273099467, 174110.85507744317874312 179559.41203175002010539))"
    thematic_dict = {"id_1": geom_from_wkt(wkt)}
    loader = DictLoader(data_dict=thematic_dict)

    aligner.load_thematic_data(loader)
    # Load reference data
    loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
    aligner.load_reference_data(loader)

    dict_evaluated2, prop_dictionary2 = aligner.evaluate()

    wkt = "LINESTRING (171741.11190000033820979 171839.01070547936251387, 171751.68948904142598622 171847.23771917796693742, 171762.26707808251376264 171855.69979041084297933, 171772.72713835645117797 171862.9865739724773448, 171783.53978493178146891 171870.8610013697471004, 171794.94007534274714999 171877.79519862998859026, 171801.87427260301774368 171884.2592808217741549)"
    thematic_dict = {"id_1": geom_from_wkt(wkt)}
    loader = DictLoader(data_dict=thematic_dict)

    aligner.load_thematic_data(loader)
    # Load reference data
    loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
    aligner.load_reference_data(loader)

    dict_evaluated3, prop_dictionary3 = aligner.evaluate()

    print(dict_evaluated1 == dict_evaluated3)
    print("END LINSTRING")
