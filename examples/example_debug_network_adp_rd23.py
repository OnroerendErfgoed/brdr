from shapely import from_wkt

from brdr.aligner import Aligner
from brdr.be.grb.enums import GRBType
from brdr.be.grb.loader import GRBActualLoader
from brdr.configs import AlignerConfig, ProcessorConfig
from brdr.constants import REMARK_FIELD_NAME
from brdr.loader import DictLoader
from brdr.processor import NetworkGeometryProcessor


if __name__ == "__main__":
    # Input polygon from issue case
    polygon_wkt = (
        "POLYGON ((173059.31170800820109434 171272.14827762849745341, "
        "173057.21826837002299726 171280.74744415286113508, "
        "173057.14863518998026848 171281.03327300003729761, "
        "173067.43548318991088308 171283.63436100000399165, "
        "173069.67049120008596219 171284.19948100007604808, "
        "173071.05097119999118149 171279.19135298999026418, "
        "173070.42377120009041391 171279.02226498996606097, "
        "173071.18370720007806085 171275.86220098999910988, "
        "173067.18370719999074936 171274.69228098992607556, "
        "173067.33973920001881197 171274.10220098998979665, "
        "173064.79074719990603626 171273.48127299008774571, "
        "173065.64194718989892863 171269.9832889900135342, "
        "173060.16303518999484368 171268.6503609899955336, "
        "173059.31170800820109434 171272.14827762849745341))"
    )
    thematic_geom = from_wkt(polygon_wkt)

    processor_config = ProcessorConfig(
        od_strategy=ProcessorConfig().od_strategy,
        network_use_directed_graph=False,
    )
    processor = NetworkGeometryProcessor(config=processor_config)

    aligner = Aligner(
        crs="EPSG:31370",
        processor=processor,
        config=AlignerConfig(profile_performance=True),
        feedback=None,
    )

    aligner.load_thematic_data(DictLoader({"debug_poly": thematic_geom}))
    aligner.load_reference_data(
        GRBActualLoader(
            grb_type=GRBType.ADP,
            partition=1000,
            aligner=aligner,
        )
    )

    rd = 23.0
    result = aligner.process(relevant_distances=[rd]).results["debug_poly"][rd]

    out_geom = result["result"]
    diff = result["result_diff"]
    remarks = result.get("properties", {}).get(REMARK_FIELD_NAME, [])

    print("\n=== DEBUG SUMMARY ===")
    print(f"relevant_distance = {rd}")
    print(f"input area        = {thematic_geom.area:.3f}")
    print(f"result area       = {out_geom.area:.3f}")
    print(f"diff area         = {diff.area:.3f}")
    print(f"remarks           = {remarks}")
    print(f"input centroid    = {thematic_geom.centroid.wkt}")
    print(f"result centroid   = {out_geom.centroid.wkt}")
    print("\nInput WKT:")
    print(thematic_geom.wkt)
    print("\nResult WKT:")
    print(out_geom.wkt)
