from pyproj import Transformer
from shapely.ops import transform

from brdr.aligner import Aligner
from brdr.be.grb.enums import GRBType
from brdr.be.grb.loader import GRBActualLoader
from brdr.be.oe.enums import OEType
from brdr.be.oe.loader import OnroerendErfgoedLoader
from brdr.loader import DictLoader
from brdr.viz import show_map


def transforms_geom_31370_to_3812(geom):
    """
    Transforms a Shapely (Multi-)Polygon from EPSG:31370 to EPSG:3812.
    """
    transformer = Transformer.from_crs("EPSG:31370", "EPSG:3812", always_xy=True)

    return transform(transformer.transform, geom)


if __name__ == "__main__":
    # 1) Load data in Lambert72  & Align to reference layer (72)
    ##################################################################

    # BASE
    # =====
    # Initiate an Aligner to create a dataset that is base-referenced on a specific base_year
    aligner72 = Aligner(crs="EPSG:31370")
    print("start loading OE-objects")
    # Load the thematic data to evaluate
    loader = OnroerendErfgoedLoader(
        objectids=["https://id.erfgoed.net/aanduidingsobjecten/121125"],
        oetype=OEType.AO,
    )
    aligner72.load_thematic_data(loader)
    print(
        "Number of OE-thematic features loaded into base-aligner: "
        + str(len(aligner72.thematic_data.features))
    )
    aligner72.load_reference_data(
        GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner72)
    )
    print("Reference-data loaded")
    #
    rd = 2
    aligner_result72 = aligner72.process(relevant_distances=[rd])
    print(
        "Processed: the data-object is now fully aligned to the reference-data in Lambert72"
    )
    dict08 = {}

    # transform features (72-->2008)
    print(
        "Next, we do a transformation of the data-object to Lambert08, possibly resulting in some deviations"
    )
    for key, processresult in aligner_result72.results.items():
        dict08[key] = transforms_geom_31370_to_3812(processresult[rd]["result"])

    # Align to GRB (2008)
    print(
        "We make a new aligner with CRS EPSG:3812 (Lambert08), so we can align the deviated object to the reference-data(Lambert08)"
    )
    aligner08 = Aligner(crs="EPSG:3812")
    loader = DictLoader(dict08)
    aligner08.load_thematic_data(loader)
    aligner08.load_reference_data(
        GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner08)
    )
    rd08 = 0.5  # we use 0.5, so deviations up to 50cm will be aligner to the reference-data(Lambert08)
    aligner_result08 = aligner08.process(relevant_distances=[rd08])
    for key, processresult in aligner_result08.results.items():

        print(
            "Resulting aligned geometry (Lambert08): "
            + processresult[rd08]["result"].wkt
        )
        print(
            "Corrected deviations (Lambert08): "
            + processresult[rd08]["result_diff"].wkt
        )

    thematic_geometries = {
        key: feat.geometry for key, feat in aligner08.thematic_data.features.items()
    }
    reference_geometries = {
        key: feat.geometry for key, feat in aligner08.reference_data.features.items()
    }
    show_map(aligner_result08.results, thematic_geometries, reference_geometries)
