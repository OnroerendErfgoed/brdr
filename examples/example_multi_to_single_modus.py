from brdr.aligner import Aligner
from brdr.be.grb.enums import GRBType
from brdr.be.grb.loader import GRBActualLoader
from brdr.be.oe.enums import OEType
from brdr.be.oe.loader import OnroerendErfgoedLoader
from brdr.configs import AlignerConfig
from brdr.enums import OpenDomainStrategy
from brdr.viz import print_observation_of_aligner_results
from brdr.viz import show_map

if __name__ == "__main__":
    """
    # This example shows the usage of the setting 'multi_as_single_modus'
    # True: (default): All polygons inside a MultiPolygon will be processed seperataly by the algorithm, and merged after processing.
    # False: Multipolygon will be processed directly by the algorithm

    """
    # Example (ErfgoedObject): https://inventaris.onroerenderfgoed.be/erfgoedobjecten/305858
    loader = OnroerendErfgoedLoader([305858], oetype=OEType.EO)
    relevant_distance = 5  # rd is taken very high to show the difference
    od_strategy = OpenDomainStrategy.SNAP_ALL_SIDE
    threshold_circle_ratio = 0.75  # default it is 0.98, but because it are not fully circles in this example we put this on 0.75

    # EXAMPLE of "multi_as_single_modus"=FALSE
    print("EXAMPLE with 'multi_as_single_modus'=False")
    aligner_config = AlignerConfig()
    aligner_config.multi_as_single_modus = False
    aligner = Aligner(
        config=aligner_config,
    )
    aligner.load_thematic_data(loader)
    aligner.load_reference_data(
        GRBActualLoader(aligner=aligner, grb_type=GRBType.ADP, partition=1000)
    )
    aligner_result = aligner.process(relevant_distances=[relevant_distance])
    aligner_result.save_results(path="output/", aligner=aligner)
    print_observation_of_aligner_results(
        aligner_result.get_results(aligner=aligner), aligner
    )
    thematic_geometries = {
        key: feat.geometry for key, feat in aligner.thematic_data.features.items()
    }
    reference_geometries = {
        key: feat.geometry for key, feat in aligner.reference_data.features.items()
    }
    show_map(aligner_result.results, thematic_geometries, reference_geometries)

    # WITH "multi_as_single_modus"=True
    print("EXAMPLE with 'multi_as_single_modus'=True")
    aligner_config = AlignerConfig()
    aligner_config.multi_as_single_modus = True
    aligner = Aligner(
        config=aligner_config,
    )
    aligner.load_thematic_data(loader)
    aligner.load_reference_data(
        GRBActualLoader(aligner=aligner, grb_type=GRBType.ADP, partition=1000)
    )
    aligner_result = aligner.process(relevant_distances=[relevant_distance])
    aligner_result.save_results(path="output/", aligner=aligner)

    print_observation_of_aligner_results(
        aligner_result.get_results(aligner=aligner), aligner
    )

    thematic_geometries = {
        key: feat.geometry for key, feat in aligner.thematic_data.features.items()
    }
    reference_geometries = {
        key: feat.geometry for key, feat in aligner.reference_data.features.items()
    }
    show_map(aligner_result.results, thematic_geometries, reference_geometries)
