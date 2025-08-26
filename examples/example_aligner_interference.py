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

    wkt = "Polygon ((174906.57643806317355484 179830.59888437716290355, 174719.95761857370962389 179820.51138062097015791, 174538.38255096235661767 179691.89570772959268652, 174442.55126527859829366 179555.71440702106337994, 174364.37311116815544665 179432.14248600779683329, 174576.21069004805758595 179376.66121534878038801, 174783.00451704987790436 179394.31434692209586501, 174906.57643806317355484 179830.59888437716290355))"
    thematic_dict = {"id_1": geom_from_wkt(wkt)}
    loader = DictLoader(data_dict=thematic_dict)

    aligner.load_thematic_data(loader)
    # Load reference data
    loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
    aligner.load_reference_data(loader)

    # PREDICT the 'stable' relevant distances, for a series of relevant distances

    # predict which relevant distances are interesting to propose as resulting geometry
    dict_evaluated1, prop_dictionary1 = aligner.evaluate()

    wkt = "Polygon ((174015.08694170592934825 179025.39916784031083807, 174040.71934808720834553 179031.93985084796440788, 174037.1838437587430235 178986.50862022733781487, 174030.64316075111855753 178982.97311589887249283, 174016.85469387014745735 179002.24161448894301429, 174021.62762471355381422 179005.77711881740833633, 174018.62244603436556645 179008.60552228015149012, 174015.08694170592934825 179025.39916784031083807))"
    thematic_dict = {"id_1": geom_from_wkt(wkt)}
    loader = DictLoader(data_dict=thematic_dict)

    aligner.load_thematic_data(loader)
    # Load reference data
    loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
    aligner.load_reference_data(loader)

    dict_evaluated2, prop_dictionary2 = aligner.evaluate()

    wkt = "Polygon ((174906.57643806317355484 179830.59888437716290355, 174719.95761857370962389 179820.51138062097015791, 174538.38255096235661767 179691.89570772959268652, 174442.55126527859829366 179555.71440702106337994, 174364.37311116815544665 179432.14248600779683329, 174576.21069004805758595 179376.66121534878038801, 174783.00451704987790436 179394.31434692209586501, 174906.57643806317355484 179830.59888437716290355))"
    thematic_dict = {"id_1": geom_from_wkt(wkt)}
    loader = DictLoader(data_dict=thematic_dict)

    aligner.load_thematic_data(loader)
    # Load reference data
    loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
    aligner.load_reference_data(loader)

    dict_evaluated3, prop_dictionary3 = aligner.evaluate()

    print(dict_evaluated1 == dict_evaluated3)
    print("END POLYGON")

    "LINESTRING (171741.11190000033820979 171839.01070547936251387, 171751.68948904142598622 171847.23771917796693742, 171762.26707808251376264 171855.69979041084297933, 171772.72713835645117797 171862.9865739724773448, 171783.53978493178146891 171870.8610013697471004, 171794.94007534274714999 171877.79519862998859026, 171801.87427260301774368 171884.2592808217741549)"
