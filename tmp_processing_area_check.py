from shapely import from_wkt
from brdr.aligner import Aligner
from brdr.loader import DictLoader
from brdr.geometry_utils import safe_symmetric_difference


def check_case(name, thematic_wkt, reference_wkt, scope_wkt):
    a = Aligner()
    thematic = {'t1': from_wkt(thematic_wkt)}
    reference = {'r1': from_wkt(reference_wkt)}
    scope = from_wkt(scope_wkt)
    a.load_thematic_data(DictLoader(thematic))
    a.load_reference_data(DictLoader(reference))
    pr = a.process([1], processing_area=scope).results['t1'][1]
    result = pr.get('result')
    diff = pr.get('result_diff')
    orig = thematic['t1']
    expected = safe_symmetric_difference(result, orig)
    ok_geom = result is not None
    if diff is None:
        ok_diff = expected is None or expected.is_empty
        d1 = 0.0
    else:
        d1 = getattr(diff, 'area', 0.0)
        d2 = 0.0 if expected is None else getattr(expected, 'area', 0.0)
        ok_diff = round(abs(d1-d2),6) == 0
    return name, ok_geom, ok_diff

cases = [
    ('partial_scope_polygon','POLYGON ((0 0, 0 6, 6 6, 6 0, 0 0))','POLYGON ((0 1, 0 7, 7 7, 7 1, 0 1))','POLYGON ((0 0, 0 6, 3 6, 3 0, 0 0))'),
    ('partial_scope_line','LINESTRING (0 0, 10 0)','LINESTRING (0 1, 10 1)','POLYGON ((0 -1, 0 1, 5 1, 5 -1, 0 -1))'),
    ('no_overlap','LINESTRING (0 0, 10 0)','LINESTRING (0 1, 10 1)','POLYGON ((20 20, 20 21, 21 21, 21 20, 20 20))')
]

all_ok = True
for c in cases:
    name, ok_geom, ok_diff = check_case(*c)
    ok = ok_geom and ok_diff
    all_ok = all_ok and ok
    print(f"{name}: ok={ok} ok_geom={ok_geom} ok_diff={ok_diff}")

raise SystemExit(0 if all_ok else 1)
