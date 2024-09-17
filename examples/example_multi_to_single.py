from brdr.aligner import Aligner
from brdr.enums import GRBType
from brdr.grb import GRBActualLoader
from brdr.oe import OnroerendErfgoedLoader
from examples import print_formula
from examples import show_map

# EXAMPLE of "multi_as_single_modus"

# Initiate brdr
aligner = Aligner()
# WITHOUT MULTI_TO_SINGLE
aligner.multi_as_single_modus = False
# Load thematic data & reference data
# Get a specific feature of OE that exists out of a Multipolygon
loader = OnroerendErfgoedLoader([110082])
aligner.load_thematic_data(loader)
aligner.load_reference_data(
    GRBActualLoader(aligner=aligner, grb_type=GRBType.GBG, partition=1000)
)

rel_dist = 20
dict_results_by_distance = {rel_dist: aligner.process_dict_thematic(rel_dist, 4)}
aligner.export_results("output/")
show_map(dict_results_by_distance, aligner.dict_thematic, aligner.dict_reference)

print_formula(dict_results_by_distance, aligner)

# WITH MULTI_TO_SINGLE

# Initiate brdr
aligner = Aligner()
aligner.multi_as_single_modus = True
# Load thematic data & reference data
# Get a specific feature of OE that exists out of a Multipolygon
loader = OnroerendErfgoedLoader([110082])
aligner.load_thematic_data(loader)
aligner.load_reference_data(
    GRBActualLoader(aligner=aligner, grb_type=GRBType.GBG, partition=1000)
)

rel_dist = 20
dict_results_by_distance = {rel_dist: aligner.process_dict_thematic(rel_dist, 4)}
aligner.export_results("output/")
show_map(dict_results_by_distance, aligner.dict_thematic, aligner.dict_reference)

print_formula(dict_results_by_distance, aligner)
