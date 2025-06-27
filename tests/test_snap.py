import unittest

from shapely import GeometryCollection

from brdr.enums import SnapStrategy
from brdr.geometry_utils import (
    geom_from_wkt,
    snap_geometry_to_reference,
    safe_unary_union,
    geometric_equality,
)


class TestSnap(unittest.TestCase):
    """
    Testcases for snapping geometries based on snap_geometry_to_reference.
    Testdata is defined here
    """

    point_thematic = geom_from_wkt(
        "POINT (174043.75093678693519905 179298.98393741055042483)"
    )
    point_reference = geom_from_wkt(
        "POINT (174079.35689694649772719 179283.04217797549790703)"
    )
    line_thematic = geom_from_wkt(
        "LINESTRING (174021.41628905048128217 179291.53905483175185509, 174081.29904022792470641 179294.45226975387777202)"
    )
    line_reference = geom_from_wkt(
        "LINESTRING (174039.13834649353520945 179282.79941006531589665, 174038.16727485283627175 179272.92684838469722308, 174052.65242682682583109 179273.33146156833390705, 174093.43743573685060255 179285.46985707725980319)"
    )
    polygon_thematic = geom_from_wkt(
        "POLYGON ((174029.67039799655321985 179306.67158789956010878, 174048.44444971703342162 179306.18605207919608802, 174048.92998553739744239 179258.92723223107168451, 174030.15593381691724062 179259.41276805143570527, 174029.67039799655321985 179306.67158789956010878))"
    )
    polygon_reference = geom_from_wkt(
        "POLYGON ((174054.91826065513305366 179305.53867098540649749, 174056.37486811622511595 179258.76538695761701092, 174072.72124073491431773 179260.38383969213464297, 174070.61725218003266491 179306.99527844646945596, 174054.91826065513305366 179305.53867098540649749))"
    )

    point_t1 = geom_from_wkt(
        "POINT (174043.75093678693519905 179298.98393741055042483)"
    )
    point_t2 = geom_from_wkt(
        "POINT (174079.35689694649772719 179283.04217797549790703)"
    )
    line_t1 = geom_from_wkt(
        "LINESTRING (174021.41628905048128217 179291.53905483175185509, 174081.29904022792470641 179294.45226975387777202)"
    )
    line_t2 = geom_from_wkt(
        "LINESTRING (174039.13834649353520945 179282.79941006531589665, 174038.16727485283627175 179272.92684838469722308, 174052.65242682682583109 179273.33146156833390705, 174093.43743573685060255 179285.46985707725980319)"
    )
    polygon_t1 = geom_from_wkt(
        "POLYGON ((174029.67039799655321985 179306.67158789956010878, 174048.44444971703342162 179306.18605207919608802, 174048.92998553739744239 179258.92723223107168451, 174030.15593381691724062 179259.41276805143570527, 174029.67039799655321985 179306.67158789956010878))"
    )
    polygon_t2 = geom_from_wkt(
        "POLYGON ((174054.91826065513305366 179305.53867098540649749, 174056.37486811622511595 179258.76538695761701092, 174072.72124073491431773 179260.38383969213464297, 174070.61725218003266491 179306.99527844646945596, 174054.91826065513305366 179305.53867098540649749))"
    )
    geometries_thematic = [point_t1, point_t2, line_t1, line_t2, polygon_t1, polygon_t2]
    geometrycollection_thematic = GeometryCollection(geometries_thematic)

    point_r1 = geom_from_wkt(
        "POINT (174043.75093678693519905 179298.98393741055042483)"
    )
    point_r2 = geom_from_wkt(
        "POINT (174079.35689694649772719 179283.04217797549790703)"
    )
    line_r1 = geom_from_wkt(
        "LINESTRING (174021.41628905048128217 179291.53905483175185509, 174081.29904022792470641 179294.45226975387777202)"
    )
    line_r2 = geom_from_wkt(
        "LINESTRING (174039.13834649353520945 179282.79941006531589665, 174038.16727485283627175 179272.92684838469722308, 174052.65242682682583109 179273.33146156833390705, 174093.43743573685060255 179285.46985707725980319)"
    )
    polygon_r1 = geom_from_wkt(
        "POLYGON ((174029.67039799655321985 179306.67158789956010878, 174048.44444971703342162 179306.18605207919608802, 174048.92998553739744239 179258.92723223107168451, 174030.15593381691724062 179259.41276805143570527, 174029.67039799655321985 179306.67158789956010878))"
    )
    polygon_r2 = geom_from_wkt(
        "POLYGON ((174054.91826065513305366 179305.53867098540649749, 174056.37486811622511595 179258.76538695761701092, 174072.72124073491431773 179260.38383969213464297, 174070.61725218003266491 179306.99527844646945596, 174054.91826065513305366 179305.53867098540649749))"
    )
    geometries_reference = [
        point_r1,
        point_r2,
        line_r1,
        line_r2,
        polygon_r1,
        polygon_r2,
    ]

    geometrycollection_reference = GeometryCollection(geometries_reference)

    def test_snap_point_to_polygon(self):
        """
        Test to snap point to polygon reference
        :return:
        """
        point_thematic = self.point_thematic
        geom_reference = self.polygon_reference
        geom_snapped = snap_geometry_to_reference(
            point_thematic,
            geom_reference,
            snap_strategy=SnapStrategy.PREFER_VERTICES,
            max_segment_length=2,
            tolerance=5,
        )
        assert geom_snapped.geom_type == point_thematic.geom_type

    def test_snap_line_to_polygon(self):
        """
        Test to snap line to polygon reference
        :return:
        """
        line_thematic = self.line_thematic
        geom_reference = self.polygon_reference
        geom_snapped = snap_geometry_to_reference(
            line_thematic,
            geom_reference,
            snap_strategy=SnapStrategy.PREFER_VERTICES,
            max_segment_length=2,
            tolerance=5,
        )
        assert geom_snapped.geom_type == line_thematic.geom_type

    def test_snap_polygon_to_polygon(self):
        """
        Test to snap polygon to polygon reference
        :return:
        """
        polygon_thematic = self.polygon_thematic
        geom_reference = self.polygon_reference
        geom_snapped = snap_geometry_to_reference(
            polygon_thematic,
            geom_reference,
            snap_strategy=SnapStrategy.PREFER_VERTICES,
            max_segment_length=2,
            tolerance=5,
        )
        assert geom_snapped.geom_type == polygon_thematic.geom_type

    def test_snap_geometrycollection_to_polygon(self):
        """
        Test to snap geometrycollection to polygon reference
        :return:
        """
        geometrycollection_thematic = self.geometrycollection_thematic
        geom_reference = self.polygon_reference
        geom_snapped = snap_geometry_to_reference(
            geometrycollection_thematic,
            geom_reference,
            snap_strategy=SnapStrategy.PREFER_VERTICES,
            max_segment_length=2,
            tolerance=5,
        )
        assert geom_snapped.geom_type == geometrycollection_thematic.geom_type

    def test_snap_point_to_line(self):
        """
        Test to snap point to line reference
        :return:
        """
        point_thematic = self.point_thematic
        geom_reference = safe_unary_union([self.line_thematic, self.line_reference])

        geom_snapped = snap_geometry_to_reference(
            point_thematic,
            geom_reference,
            snap_strategy=SnapStrategy.PREFER_VERTICES,
            max_segment_length=2,
            tolerance=5,
        )
        assert geom_snapped.geom_type == point_thematic.geom_type

    def test_snap_line_to_line(self):
        """
        Test to snap line to line reference
        :return:
        """
        line_1 = self.line_thematic
        geom_reference = safe_unary_union([self.line_thematic, self.line_reference])

        geom_snapped = snap_geometry_to_reference(
            line_1,
            geom_reference,
            snap_strategy=SnapStrategy.PREFER_VERTICES,
            max_segment_length=2,
            tolerance=5,
        )
        assert geom_snapped.geom_type == line_1.geom_type

    def test_snap_polygon_to_line(self):
        """
        Test to snap polygon to polygon reference
        :return:
        """
        polygon_1 = self.polygon_thematic
        geom_reference = safe_unary_union([self.line_thematic, self.line_reference])

        geom_snapped = snap_geometry_to_reference(
            polygon_1,
            geom_reference,
            snap_strategy=SnapStrategy.PREFER_VERTICES,
            max_segment_length=2,
            tolerance=5,
        )
        assert geom_snapped.geom_type == polygon_1.geom_type

    def test_snap_geometrycollection_to_line(self):
        """
        Test to snap geometrycollection to line reference
        :return:
        """
        geom_thematic = self.geometrycollection_thematic
        geom_reference = self.line_reference
        geom_snapped = snap_geometry_to_reference(
            geom_thematic,
            geom_reference,
            snap_strategy=SnapStrategy.PREFER_VERTICES,
            max_segment_length=2,
            tolerance=5,
        )
        assert geom_snapped.geom_type == geom_thematic.geom_type

    ##############
    def test_snap_point_to_point(self):
        """
        Test to snap point to point reference
        :return:
        """
        point_1 = self.point_thematic
        geom_reference = self.point_reference

        geom_snapped = snap_geometry_to_reference(
            point_1,
            geom_reference,
            snap_strategy=SnapStrategy.PREFER_VERTICES,
            max_segment_length=2,
            tolerance=5,
        )
        assert geom_snapped.geom_type == point_1.geom_type

    def test_snap_line_to_point(self):
        """
        Test to snap line to point reference
        :return:
        """
        line_1 = self.line_thematic
        geom_reference = self.point_reference
        geom_snapped = snap_geometry_to_reference(
            line_1,
            geom_reference,
            snap_strategy=SnapStrategy.PREFER_VERTICES,
            max_segment_length=2,
            tolerance=5,
        )
        assert geom_snapped.geom_type == line_1.geom_type

    def test_snap_polygon_to_point(self):
        """
        Test to snap polygon to polygon reference
        :return:
        """
        polygon_1 = self.polygon_thematic
        geom_reference = self.point_reference
        geom_snapped = snap_geometry_to_reference(
            polygon_1,
            geom_reference,
            snap_strategy=SnapStrategy.PREFER_VERTICES,
            max_segment_length=2,
            tolerance=40,
        )
        assert geom_snapped.geom_type == polygon_1.geom_type

    def test_snap_polygon_to_point_no_preference(self):
        """
        Test to snap polygon to polygon reference
        :return:
        """
        polygon_1 = geom_from_wkt(
            "MULTIPOLYGON (((171795.71618631482124329 171817.88460136577486992, 171784.53532230854034424 171806.1688893586397171, 171746.73993028700351715 171841.20300138369202614, 171746.28380228579044342 171841.62578538432717323, 171767.35881029814481735 171856.89906539395451546, 171767.47471430152654648 171856.76376939192414284, 171798.38581831753253937 171820.68191336840391159, 171798.01820231974124908 171820.29669736698269844, 171795.71618631482124329 171817.88460136577486992)))"
        )
        geom_reference = geom_from_wkt(
            "MULTIPOINT ((171756.52506366037414409 171850.34457502837176435),(171766.12778689494007267 171857.04121096828021109),(171777.13617191879893653 171865.19089055553195067),(171745.72200002145837061 171841.94219219812657684),(171770.11351679253857583 171840.12903735807049088))"
        )
        geom_snapped = snap_geometry_to_reference(
            polygon_1,
            geom_reference,
            snap_strategy=SnapStrategy.NO_PREFERENCE,
            max_segment_length=2,
            tolerance=5,
        )
        assert not geometric_equality(polygon_1, geom_snapped, 0.01, mitre_limit=5)

    def test_snap_geometrycollection_to_point(self):
        """
        Test to snap geometrycollection to point reference
        :return:
        """
        geom_thematic = self.geometrycollection_thematic
        geom_reference = geom_from_wkt(
            "POINT (174028.67039799655321985 179306.67158789956010878)"
        )
        geom_snapped = snap_geometry_to_reference(
            geom_thematic,
            geom_reference,
            snap_strategy=SnapStrategy.PREFER_VERTICES,
            max_segment_length=2,
            tolerance=5,
        )
        assert geom_snapped.geom_type == geom_thematic.geom_type

    def test_snap_point_to_geometrycollection(self):
        """
        Test to snap point to geometrycollection reference
        :return:
        """
        point_1 = self.point_thematic
        geom_reference = self.geometrycollection_reference
        geom_snapped = snap_geometry_to_reference(
            point_1,
            geom_reference,
            snap_strategy=SnapStrategy.PREFER_VERTICES,
            max_segment_length=2,
            tolerance=5,
        )
        assert geom_snapped.geom_type == point_1.geom_type

    def test_snap_line_to_geometrycollection(self):
        """
        Test to snap line to geometrycollection reference
        :return:
        """
        line_1 = self.line_thematic
        geom_reference = self.geometrycollection_reference
        geom_snapped = snap_geometry_to_reference(
            line_1,
            geom_reference,
            snap_strategy=SnapStrategy.PREFER_VERTICES,
            max_segment_length=2,
            tolerance=5,
        )

        assert geom_snapped.geom_type == line_1.geom_type

    def test_snap_polygon_to_geometrycollection(self):
        """
        Test to snap polygon to geometrycollection reference
        :return:
        """
        polygon_1 = self.polygon_thematic
        geom_reference = self.geometrycollection_reference
        geom_snapped = snap_geometry_to_reference(
            polygon_1,
            geom_reference,
            snap_strategy=SnapStrategy.PREFER_VERTICES,
            max_segment_length=2,
            tolerance=5,
        )
        assert geom_snapped.geom_type == polygon_1.geom_type

    def test_snap_geometrycollection_to_geometrycollection(self):
        """
        Test to snap geometrycollection to geometrycollection reference
        :return:
        """
        geom_thematic = self.geometrycollection_thematic
        geom_reference = self.geometrycollection_reference
        geom_snapped = snap_geometry_to_reference(
            geom_thematic,
            geom_reference,
            snap_strategy=SnapStrategy.PREFER_VERTICES,
            max_segment_length=2,
            tolerance=5,
        )
        assert geom_snapped.geom_type == geom_thematic.geom_type
