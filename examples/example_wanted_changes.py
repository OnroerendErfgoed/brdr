import numpy as np
from brdr.aligner import Aligner
from brdr.utils import get_breakpoints
from examples import show_results, plot_series, show_individual_results

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    """
    example to test if we can seperate wanted deviations from unwanted deviations.
    By calculating results with a series of relevant distances we get a graphic.
    There is a difference in the graphic for wanted and unwanted deviations:
    * unwanted deviations: will evoluate more gradually
    * wanted deviations: have a harder breakpoint
    The type of breakpoint, breakpoint-distance and length between breakpoints can gave an indication of wanted and unwanted deviations
    """
    # TODO: future possibilities to detect the 'wanted' breakpoint automatically? (AI-machine learning?)
    ##Initiate a Aligner
    aligner = Aligner()
    ##Load thematic data & reference data
    aligner.load_thematic_data_file("../tests/testdata/test_wanted_changes.geojson", 'theme_id')
    aligner.load_reference_data_grb_actual(grb_type='adp', partition=1000) #gebruik de actuele adp-percelen adp= administratieve percelen



    #Example how to use the Aligner
    r,rd, rd_plus,rd_min,sd,si = aligner.process_dict_thematic(2, 4)
    aligner.export_results("output/")
    show_results(r,rd_plus,rd_min, aligner.dict_thematic, aligner.dict_reference)

    #Possibility to get the descriptive formula of a thematic feature
    for key in r:
        aligner.get_formula(r[key])

    #Example how to use a series (for histogram)
    #series = [0.2, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 6, 8, 10]
    #series = [0.1,0.2,0.3,0.4, 0.5,1,2,3,4]
    series = np.arange(0.1, 10.05, 0.1, dtype=float)
    resulting_areas = aligner.process_series(series,4,50)
    plot_series(series, resulting_areas)
    for key in resulting_areas:
        if len(resulting_areas[key]) == len(series):
            lst_diffs = list(resulting_areas[key].values())
            extremen, zero_streak = get_breakpoints(series, lst_diffs)
            print (str(key))
            for extremum in extremen:
                print(f"{extremum[0]:.2f}, {extremum[1]:.2f} ({extremum[2]})")
            for st in zero_streak:
                print(f"{st[0]:.2f} - {st[1]:.2f} -{st[2]:.2f} - {st[3]:.2f} - startextreme {st[4]:.2f} ")
                r, rd, rd_plus, rd_min, sd, si = aligner.process_dict_thematic(st[0], 4)
                show_individual_results(r,rd_plus, rd_min,aligner.dict_thematic,aligner.dict_reference, key)




