from shapely import from_wkt

from brdr.aligner import Aligner
from brdr.enums import GRBType
from brdr.grb import GRBActualLoader
from brdr.loader import DictLoader

if __name__ == "__main__":
    aligner = Aligner(max_workers=-1)
    thematic_dict = {
        "theme_id_1": from_wkt(
            "MULTIPOINT ((174135.22472090268274769 179531.68224664646550082),(174121.19416637552785687 179537.99109039537142962),(174121.04699272679863498 179548.5385352322482504))"
        )
    }

    loader = DictLoader(thematic_dict)
    aligner.load_thematic_data(loader)
    loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
    aligner.load_reference_data(loader)

    dict_processresults = aligner.process(relevant_distance=3)
    print(dict_processresults["theme_id_1"][3]["result"])
