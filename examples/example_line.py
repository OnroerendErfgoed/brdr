import os

from shapely import from_wkt

from brdr.aligner import Aligner
from brdr.enums import AlignerInputType, GRBType
from brdr.grb import GRBSpecificDateParcelLoader, GRBActualLoader
from brdr.loader import DictLoader
from brdr.utils import write_geojson

if __name__ == "__main__":
    aligner = Aligner(max_workers=-1)
    thematic_dict = {
        "theme_id_1": from_wkt("MULTILINESTRING ((174135.22254434687783942 179531.39341345359571278, 174121.64978590866667219 179537.66386910632718354, 174118.0780073722708039 179538.69571623904630542, 174121.25292162684490904 179549.09356042271247134, 174114.98246597408433445 179551.23662754453835078, 174116.56992310137138702 179556.95147320273099467, 174110.85507744317874312 179559.41203175002010539))")
    }

    loader = DictLoader(thematic_dict)
    aligner.load_thematic_data(loader)
    loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
    aligner.load_reference_data(loader)

    dict_processresults = aligner.process(relevant_distance=3)
    print (dict_processresults["theme_id_1"][3]['result'])