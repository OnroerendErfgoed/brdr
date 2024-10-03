from brdr.aligner import Aligner
from brdr.enums import OpenbaarDomeinStrategy, GRBType
from brdr.grb import GRBActualLoader
from brdr.loader import GeoJsonLoader
from brdr.utils import diffs_from_dict_series
from examples import plot_series
from examples import show_map

if __name__ == "__main__":
    # EXAMPLE to process a series of relevant distances
    # Initiate brdr
    aligner = Aligner()
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

    loader = GeoJsonLoader(_input=thematic_json, id_property="theme_identifier")
    aligner.load_thematic_data(loader)
    # Load reference data: The actual GRB-parcels
    aligner.load_reference_data(
        GRBActualLoader(grb_type=GRBType.ADP, partition=1000, aligner=aligner)
    )
    # PROCESS a series of relevant distances
    relevant_distances = [0.5, 1, 3, 6]
    dict_results = aligner.process(
        relevant_distances=relevant_distances,
        od_strategy=OpenbaarDomeinStrategy.SNAP_ALL_SIDE,
        threshold_overlap_percentage=50,
    )
    # SHOW results: map and plotted changes
    show_map(dict_results, aligner.dict_thematic, aligner.dict_reference)
    resulting_areas = diffs_from_dict_series(dict_results, aligner.dict_thematic)
    plot_series(relevant_distances, resulting_areas)
