from brdr.aligner import Aligner
from examples import show_results

if __name__ == "__main__":
    # Initiate brdr
    x = Aligner()

    # Load thematic data
    thematic_json = {
        "type": "FeatureCollection",
        "name": "test",
        "crs": {"type": "name", "properties": {"name": "urn:ogc:def:crs:EPSG::31370"}},
        "features": [
            {
                "type": "Feature",
                "properties": {"fid": 1100, "id": 1100, "theme_identifier": "1100"},
                "geometry": {
                    "type": "MultiPolygon",
                    "coordinates": [
                        [
                            [
                                [170020.885142877610633, 171986.324472956912359],
                                [170078.339491307124263, 172031.344329243671382],
                                [170049.567976467020344, 172070.009247593494365],
                                [170058.413533725659363, 172089.43287940043956],
                                [170071.570170061604585, 172102.403589786874363],
                                [170061.212614970601862, 172153.235667688539252],
                                [170008.670597244432429, 172137.344562214449979],
                                [169986.296310421457747, 172121.231194059771951],
                                [169979.658902874827618, 172095.676061166130239],
                                [169980.509304589155363, 172078.495172578957863],
                                [169985.424311356124235, 172065.162689057324314],
                                [169971.609697333740769, 172057.093288242496783],
                                [170020.885142877610633, 171986.324472956912359],
                            ]
                        ]
                    ],
                },
            }
        ],
    }

    x.load_thematic_data_geojson(thematic_json, "theme_identifier")
    # gebruik de actuele adp-percelen adp= administratieve percelen
    x.load_reference_data_grb_actual(grb_type="adp", partition=1000)

    # Example how to use the Aligner
    r, rd, rd_plus, rd_min, sd, si = x.process_dict_thematic(6, 4)
    out = x.get_last_version_date(x.dict_thematic["1100"])
    print(out)
    x.export_results("output/")
    show_results(r, rd)
