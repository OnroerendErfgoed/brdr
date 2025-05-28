from shapely import from_wkt

from brdr.aligner import Aligner
from brdr.enums import GRBType, OpenDomainStrategy
from brdr.grb import GRBActualLoader
from brdr.loader import DictLoader
from examples import show_geometry

wkt = "POLYGON ((172791.07652649749070406 171319.35066181665752083, 172821.75881976058008149 171308.36308382378774695, 172838.13653035368770361 171321.42378973981249146, 172879.59908881731098518 171303.80220239277696237, 172842.49009899239172228 171205.95056441865745001, 172804.55185799815808423 171219.42589591932483017, 172789.21071136664249934 171221.29171105020213872, 172757.69916693426785059 171270.63215562188997865, 172766.4063042116467841 171279.13198010693304241, 172776.56463103523128666 171296.96088024630444124, 172785.63235233756131493 171305.36084497725823894, 172782.68072999455034733 171312.79389392648590729, 172791.07652649749070406 171319.35066181665752083))"
theme_id = "theme_id_1"
thematic_dict = {theme_id: from_wkt(wkt)}
sample_aligner=Aligner()
sample_aligner.load_thematic_data(DictLoader(thematic_dict))
# LOAD REFERENCE DICTIONARY
sample_aligner.load_reference_data(
    GRBActualLoader(
        aligner=sample_aligner, grb_type=GRBType.ADP, partition=1000
    )
)
relevant_distance = 1
result_dict = sample_aligner.process(
    od_strategy=OpenDomainStrategy.SNAP_ALL_SIDE, relevant_distance=relevant_distance
)
show_geometry(result_dict[theme_id][relevant_distance].get("result"))
print("done")
