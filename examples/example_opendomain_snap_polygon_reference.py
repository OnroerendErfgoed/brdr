from brdr.aligner import Aligner
from brdr.configs import ProcessorConfig
from brdr.enums import OpenDomainStrategy
from brdr.loader import DictLoader
from brdr.processor import SnapGeometryProcessor
from shapely import from_wkt


if __name__ == "__main__":
    """Example: compare OD strategies in SnapGeometryProcessor with polygon references."""

    thematic = {
        "poly_1": from_wkt("POLYGON ((0 0, 0 6, 6 6, 6 0, 0 0))"),
    }
    reference = {
        "ref_poly_1": from_wkt("POLYGON ((1 1, 1 5, 5 5, 5 1, 1 1))"),
    }
    relevant_distance = 0.5
    strategies = [
        OpenDomainStrategy.EXCLUDE,
        OpenDomainStrategy.AS_IS,
        OpenDomainStrategy.SNAP_ALL_SIDE,
    ]

    for strategy in strategies:
        processor = SnapGeometryProcessor(config=ProcessorConfig(od_strategy=strategy))
        aligner = Aligner(processor=processor)
        aligner.load_thematic_data(DictLoader(thematic))
        aligner.load_reference_data(DictLoader(reference))

        result = aligner.process([relevant_distance]).results["poly_1"][relevant_distance][
            "result"
        ]
        print(
            f"{strategy.value}: type={result.geom_type}, empty={result.is_empty}, area={round(result.area, 3)}"
        )
