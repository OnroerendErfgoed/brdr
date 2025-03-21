from brdr.aligner import Aligner
from brdr.enums import GRBType, OpenDomainStrategy
from brdr.grb import GRBActualLoader
from brdr.oe import OnroerendErfgoedLoader, OEType
from examples import print_brdr_formula
from examples import show_map

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
    aligner = Aligner(
        multi_as_single_modus=False,
        relevant_distance=relevant_distance,
        od_strategy=od_strategy,
        threshold_circle_ratio=threshold_circle_ratio,
    )
    aligner.load_thematic_data(loader)
    aligner.load_reference_data(
        GRBActualLoader(aligner=aligner, grb_type=GRBType.ADP, partition=1000)
    )
    dict_results = aligner.process()
    aligner.save_results("output/")
    print_brdr_formula(dict_results, aligner)
    show_map(dict_results, aligner.dict_thematic, aligner.dict_reference)

    # WITH "multi_as_single_modus"=True
    print("EXAMPLE with 'multi_as_single_modus'=True")
    aligner = Aligner(
        multi_as_single_modus=True,
        relevant_distance=relevant_distance,
        od_strategy=od_strategy,
        threshold_circle_ratio=threshold_circle_ratio,
    )
    aligner.load_thematic_data(loader)
    aligner.load_reference_data(
        GRBActualLoader(aligner=aligner, grb_type=GRBType.ADP, partition=1000)
    )
    dict_results = aligner.process()
    aligner.save_results("output/")
    print_brdr_formula(dict_results, aligner)
    show_map(dict_results, aligner.dict_thematic, aligner.dict_reference)
