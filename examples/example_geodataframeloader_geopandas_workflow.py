import geopandas as gpd
from shapely import from_wkt

from brdr.aligner import Aligner
from brdr.loader import GeoDataFrameLoader


if __name__ == "__main__":
    """Example: use GeoDataFrameLoader for thematic/reference and export results as GeoDataFrame."""

    thematic_gdf = gpd.GeoDataFrame(
        {
            "theme_id": ["t1", "t2"],
            "name": ["theme-one", "theme-two"],
            "geometry": [
                from_wkt("POLYGON ((0 0, 0 9, 5 10, 10 0, 0 0))"),
                from_wkt("POLYGON ((20 0, 20 8, 28 8, 28 0, 20 0))"),
            ],
        },
        geometry="geometry",
        crs="EPSG:31370",
    )

    reference_gdf = gpd.GeoDataFrame(
        {
            "ref_id": ["r1", "r2"],
            "source": ["synthetic", "synthetic"],
            "geometry": [
                from_wkt("POLYGON ((0 1, 0 10, 8 10, 10 1, 0 1))"),
                from_wkt("POLYGON ((21 1, 21 7, 27 7, 27 1, 21 1))"),
            ],
        },
        geometry="geometry",
        crs="EPSG:31370",
    )

    aligner = Aligner(crs="EPSG:31370")
    aligner.load_thematic_data(
        GeoDataFrameLoader(id_property="theme_id", _input=thematic_gdf)
    )
    aligner.load_reference_data(
        GeoDataFrameLoader(
            id_property="ref_id",
            _input=reference_gdf,
            is_reference=True,
        )
    )

    relevant_distances = [0.5, 1.0]
    aligner_result = aligner.process(relevant_distances=relevant_distances)

    result_gdf = aligner_result.get_results_as_geodataframe(
        aligner=aligner,
        add_metadata=False,
        add_original_attributes=True,
    )

    print("rows:", len(result_gdf))
    print("columns:", list(result_gdf.columns))
    print(
        result_gdf[
            [
                "brdr_id",
                "relevant_distance",
                "name",
                "result",
            ]
        ]
    )
