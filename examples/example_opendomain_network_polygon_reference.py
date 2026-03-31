from brdr.aligner import Aligner
from brdr.configs import ProcessorConfig
from brdr.enums import OpenDomainStrategy
from brdr.loader import DictLoader
from brdr.processor import NetworkGeometryProcessor
from shapely import from_wkt


if __name__ == "__main__":
    """Example: compare OD strategies in NetworkGeometryProcessor with polygon references."""

    thematic = {"line_1": from_wkt("LINESTRING (0 0, 10 0)")}
    reference = {"ref_poly_1": from_wkt("POLYGON ((2 -1, 8 -1, 8 1, 2 1, 2 -1))")}
    relevant_distance = 1.0

    for strategy in OpenDomainStrategy:
        processor = NetworkGeometryProcessor(
            config=ProcessorConfig(od_strategy=strategy)
        )
        aligner = Aligner(processor=processor)
        aligner.load_thematic_data(DictLoader(thematic))
        aligner.load_reference_data(DictLoader(reference))

        result = aligner.process([relevant_distance]).results["line_1"][relevant_distance][
            "result"
        ]
        print(
            f"{strategy.value}: type={result.geom_type}, empty={result.is_empty}, length={round(result.length, 3)}"
        )
