from brdr.auto_referencer import AutoReferencer
from brdr.utils import get_oe_geojson_by_ids
from brdr.utils import multipolygons_to_singles
from examples import show_results

# example to Change a dictionary form multipolygon to single before executing the
# auto_referencer. Can be used on the thematic dictionary as the reference dictionary


if __name__ == "__main__":
    # EXAMPLE for a thematic MultiPolygon
    dict_theme = get_oe_geojson_by_ids([110082])

    # WITHOUT MULTI_TO_SINGLE
    # Initiate brdr
    x = AutoReferencer()
    # Load thematic data & reference data
    # Get a specific feature of OE that exists out of a Multipolygon

    x.load_thematic_data_dict(dict_theme)
    x.load_reference_data_grb_actual(grb_type="gbg", partition=1000)

    r, rd, rd_plus, rd_min, sd, si = x.process_dict_thematic(2, 4)
    x.export_results("output/")
    show_results(r, rd)

    for key in r:
        print(key)
        print(x.get_formula(r[key]))

    # WITH MULTI_TO_SINGLE
    # Initiate brdr
    x = AutoReferencer()
    # Load thematic data & reference data
    # Get a specific feature of OE that exists out of a Multipolygon
    dict_theme = multipolygons_to_singles(dict_theme)
    x.load_thematic_data_dict(dict_theme)
    x.load_reference_data_grb_actual(grb_type="gbg", partition=1000)

    r, rd, rd_plus, rd_min, sd, si = x.process_dict_thematic(5, 4)
    x.export_results("output/")
    show_results(r, rd)

    for key in r:
        print(key)
        print(x.get_formula(r[key]))
