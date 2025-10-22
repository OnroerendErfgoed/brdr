from pyproj import Transformer
from shapely.ops import transform

from brdr.aligner import Aligner
from brdr.enums import GRBType, OpenDomainStrategy
from brdr.loaders.grb import GRBActualLoader
from brdr.loaders.loader import DictLoader
from brdr.loaders.oe import OnroerendErfgoedLoader, OEType
from examples import show_map


def transform_geom_31370_to_3812(geom):
    """
    Transformeert een Shapely Polygon of MultiPolygon van EPSG:31370 naar EPSG:3812.
    """
    transformer = Transformer.from_crs("EPSG:31370", "EPSG:3812", always_xy=True)
    return transform(transformer.transform, geom)


if __name__ == "__main__":
    # 1) Load data in Lambert72  & Align to reference layer (72)
    ##################################################################

    # BASE
    # =====
    # Initiate an Aligner to create a themeset that is base-referenced on a specific
    # base_year
    aligner72 = Aligner(crs="EPSG:31370")
    print("start loading OE-objects")
    # Load the thematic data to evaluate
    # loader = OnroerendErfgoedLoader(bbox=bbox, partition=0)
    loader = OnroerendErfgoedLoader(objectids=[121125], oetype=OEType.AO)
    aligner72.load_thematic_data(loader)
    print(
        "Number of OE-thematic features loaded into base-aligner: "
        + str(len(aligner72.dict_thematic))
    )
    aligner72.load_reference_data(
        GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner72)
    )
    print("Reference-data loaded")
    #
    rd=2
    dict_processresults = aligner72.process(relevant_distance=rd,od_strategy=OpenDomainStrategy.SNAP_ALL_SIDE)
    print("Processed")
    dict08={}
    # transform features (72-->2008)
    for key,processresult in dict_processresults.items():
        dict08[key] = transform_geom_31370_to_3812(processresult[rd]['result'])

    # Align to GRB (2008)
    aligner08= Aligner(crs="EPSG:3812")
    loader = DictLoader(dict08)
    # loader = GeoJsonLoader(id_property=id,_input=fcs)
    aligner08.load_thematic_data(loader)
    aligner08.load_reference_data(
        GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner08)
    )
    rd08=0.5
    dict_processresults08 = aligner08.process(relevant_distance=rd08,od_strategy=OpenDomainStrategy.SNAP_ALL_SIDE)
    for key,processresult in dict_processresults08.items():

        print(processresult[rd08]["result"].wkt)
        print(processresult[rd08]['result_diff'].wkt)
    show_map(dict_processresults08, aligner08.dict_thematic, aligner08.dict_reference
        )
