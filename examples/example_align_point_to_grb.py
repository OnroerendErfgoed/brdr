from shapely import from_wkt

from brdr.aligner import Aligner
from brdr.be.grb.enums import GRBType
from brdr.be.grb.loader import GRBActualLoader
from brdr.enums import AlignerResultType
from brdr.loader import DictLoader

if __name__ == "__main__":
    aligner = Aligner()
    thematic_dict = {
        "theme_id_1": from_wkt(
            "MULTIPOINT ((174135.22472090268274769 179531.68224664646550082),(174121.19416637552785687 179537.99109039537142962),(174121.04699272679863498 179548.5385352322482504))"
        )
    }

    loader = DictLoader(thematic_dict)
    aligner.load_thematic_data(loader)
    loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
    aligner.load_reference_data(loader)

    aligner_result = aligner.process(relevant_distances=[3])
    print(aligner_result.results["theme_id_1"][3]["result"])
    aligner_result = aligner.predict()
    aligner_result.save_results(
        path="output", result_type=AlignerResultType.PREDICTIONS, aligner=aligner
    )
