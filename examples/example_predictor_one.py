from brdr.aligner import Aligner
from brdr.be.grb.enums import GRBType
from brdr.be.grb.loader import GRBActualLoader
from brdr.enums import AlignerResultType, OpenDomainStrategy
from brdr.loader import GeoJsonFileLoader
from examples import show_map, plot_dict_diffs

# Press the green button in the gutter to run the script.
if __name__ == "__main__":
    """
    EXAMPLE to use the predictor-function to automatically predict which resulting
    geometries are interesting to look at (based on detection of breakpoints and
    relevant distances of 'no-change')
    """
    # Initiate an Aligner
    aligner = Aligner(max_workers=-1)
    # Load thematic data & reference data
    loader = GeoJsonFileLoader("input/predictor_one.geojson", "theme_id")

    aligner.load_thematic_data(loader)
    # Load reference data
    loader = GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
    aligner.load_reference_data(loader)

    # PREDICT the 'stable' relevant distances, for a series of relevant distances
    series = [3]
    # predict which relevant distances are interesting to propose as resulting geometry
    alignerresults = aligner.predict(
        relevant_distances=series,
    )
    dict_series = alignerresults.get_results(
    )
    # SHOW results of the predictions
    fcs_all = alignerresults.get_results_as_geojson(formula=False, aligner=aligner)

    if fcs_all is None or "result" not in fcs_all:
        print("no calculations")
    else:
        print(fcs_all["result"])
        for key in dict_series:
            plot_dict_diffs({key: aligner.diffs_dict[key]})
            show_map(
                {key: dict_series[key]},
                {key: aligner.dict_thematic[key]},
                aligner.dict_reference,
            )
