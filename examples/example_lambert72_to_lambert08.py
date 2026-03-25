from pyproj import Transformer
from shapely.ops import transform

from brdr.aligner import Aligner
from brdr.be.grb.enums import GRBType
from brdr.be.grb.loader import GRBActualLoader
from brdr.be.oe.enums import OEType
from brdr.be.oe.loader import OnroerendErfgoedLoader
from brdr.loader import DictLoader
from brdr.viz import show_map


def transform_geom_31370_to_3812(geom):
    """Transform a Shapely geometry from EPSG:31370 to EPSG:3812."""
    transformer = Transformer.from_crs("EPSG:31370", "EPSG:3812", always_xy=True)
    return transform(transformer.transform, geom)


if __name__ == "__main__":
    # 1) Load and align data in Lambert72 (EPSG:31370).
    aligner72 = Aligner(crs="EPSG:31370")
    print("start loading OE-objects")
    aligner72.load_thematic_data(
        OnroerendErfgoedLoader(
            objectids=["https://id.erfgoed.net/aanduidingsobjecten/121125"],
            oetype=OEType.AO,
        )
    )
    print(
        "Number of OE-thematic features loaded into base-aligner: "
        + str(len(aligner72.thematic_data.features))
    )
    aligner72.load_reference_data(
        GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner72)
    )
    print("Reference-data loaded")

    relevant_distance_72 = 2
    aligner_result72 = aligner72.process(relevant_distances=[relevant_distance_72])
    print(
        "Processed: the data-object is now fully aligned to the reference-data in Lambert72"
    )

    # 2) Transform aligned data from Lambert72 to Lambert08.
    dict08 = {}
    print(
        "Next, we do a transformation of the data-object to Lambert08, possibly resulting in some deviations"
    )
    for key, process_result in aligner_result72.results.items():
        dict08[key] = transform_geom_31370_to_3812(
            process_result[relevant_distance_72]["result"]
        )

    # 3) Re-align transformed data in Lambert08.
    print(
        "We make a new aligner with CRS EPSG:3812 (Lambert08), so we can align the deviated object to the reference-data(Lambert08)"
    )
    aligner08 = Aligner(crs="EPSG:3812")
    aligner08.load_thematic_data(DictLoader(dict08))
    aligner08.load_reference_data(
        GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner08)
    )

    relevant_distance_08 = 0.5
    aligner_result08 = aligner08.process(relevant_distances=[relevant_distance_08])
    for _, process_result in aligner_result08.results.items():

        print(
            "Resulting aligned geometry (Lambert08): "
            + process_result[relevant_distance_08]["result"].wkt
        )
        print(
            "Corrected deviations (Lambert08): "
            + process_result[relevant_distance_08]["result_diff"].wkt
        )

    thematic_geometries = {
        key: feat.geometry for key, feat in aligner08.thematic_data.features.items()
    }
    reference_geometries = {
        key: feat.geometry for key, feat in aligner08.reference_data.features.items()
    }
    show_map(aligner_result08.results, thematic_geometries, reference_geometries)
