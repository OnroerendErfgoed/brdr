from brdr.aligner import Aligner
from brdr.utils import get_oe_dict_by_ids
from examples import show_results

if __name__ == "__main__":
    # EXAMPLE for a thematic Polygon (aanduid_id 131635)

    # Initiate brdr
    x = Aligner()
    # Load thematic data & reference data
    dict_theme = get_oe_dict_by_ids([131635])
    x.load_thematic_data_dict(dict_theme)
    x.load_reference_data_grb_actual(grb_type="adp", partition=1000)

    #RESULTS
    r, rd, rd_plus, rd_min, sd, si = x.process_dict_thematic(2, 4)
    x.export_results("output/")
    show_results(r, rd_plus,rd_min,x.dict_thematic, x.dict_reference)
    for key in r:
        print(key)
        print(x.get_formula(r[key]))