from brdr.aligner import Aligner
from brdr.enums import GRBType
from brdr.grb import GRBActualLoader
from brdr.oe import OnroerendErfgoedLoader


# EXAMPLE to test the algorithm for erfgoedobject with relevant distance 0.2m and
# od_strategy SNAP_ALL_SIDE

# Initiate brdr
aligner = Aligner()
# Load thematic data & reference data
aanduidingsobjecten = [117798, 116800, 117881]

loader = OnroerendErfgoedLoader(aanduidingsobjecten)
aligner.load_thematic_data(loader)
loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
aligner.load_reference_data(loader)

test = aligner.process()
test = aligner.process([1, 2, 3])
test = aligner.predictor()
fcs = aligner.get_results_as_geojson(formula=True)
print(test)
print(fcs)
print(fcs["result"])
