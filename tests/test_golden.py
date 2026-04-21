import math
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
import pytest
from phylomath import bhv_geodesic, parse_newick
from phylomath.gtp import path_space_geo

def test_path_space_geo_empty():
    assert path_space_geo([], []) == 0.0

def test_path_space_geo_single_pair_is_cone_path():
    assert math.isclose(path_space_geo([3.0], [4.0]), 7.0)

def test_path_space_geo_ascending_ratios_uses_straight_line():
    assert math.isclose(path_space_geo([1.0, 3.0], [2.0, 1.0]), 5.0)

def test_path_space_geo_descending_ratios_merges_to_one_block():
    expected = math.sqrt(10) + math.sqrt(5)
    assert math.isclose(path_space_geo([3.0, 1.0], [1.0, 2.0]), expected)

def test_path_space_geo_equal_ratios_treated_as_merge_boundary():
    separate = math.sqrt((3 + 4) ** 2 + (3 + 4) ** 2)
    combined = math.sqrt(3 ** 2 + 3 ** 2) + math.sqrt(4 ** 2 + 4 ** 2)
    assert math.isclose(separate, combined)
    assert math.isclose(path_space_geo([3.0, 3.0], [4.0, 4.0]), combined)

def test_bhv_geodesic_same_tree_zero():
    t = parse_newick('((A:1,B:2):0.5,(C:3,D:0.7):0.2);')
    assert bhv_geodesic(t, t) == 0.0

def test_bhv_geodesic_same_topology_is_l2_on_shared_coords():
    t1 = parse_newick('((A:1,B:1):1,(C:1,D:1):1);')
    t2 = parse_newick('((A:2,B:2):2,(C:2,D:2):2);')
    assert math.isclose(bhv_geodesic(t1, t2), 2.0, rel_tol=1e-09)

def test_bhv_geodesic_fully_incompatible_4taxon_is_cone_path():
    t1 = parse_newick('((A:0,B:0):3,(C:0,D:0):0);')
    t2 = parse_newick('((A:0,C:0):4,(B:0,D:0):0);')
    assert math.isclose(bhv_geodesic(t1, t2), 7.0, rel_tol=1e-09)

def test_bhv_geodesic_triangle_inequality_on_5_taxa_triples():
    t1 = parse_newick('((A:1,B:1):1,(C:1,(D:1,E:1):1):1);')
    t2 = parse_newick('((A:1,C:1):1,(B:1,(D:1,E:1):1):1);')
    t3 = parse_newick('(((A:1,B:1):1,C:1):1,(D:1,E:1):1);')
    d12 = bhv_geodesic(t1, t2)
    d23 = bhv_geodesic(t2, t3)
    d13 = bhv_geodesic(t1, t3)
    assert d13 <= d12 + d23 + 1e-09
    assert d12 <= d13 + d23 + 1e-09
    assert d23 <= d12 + d13 + 1e-09

def test_bhv_geodesic_symmetry():
    t1 = parse_newick('((A:1,B:2):0.5,(C:3,D:0.7):0.2);')
    t2 = parse_newick('((A:1.1,C:0.9):0.6,(B:2.1,D:0.7):0.5);')
    assert math.isclose(bhv_geodesic(t1, t2), bhv_geodesic(t2, t1), rel_tol=1e-09)

def test_bhv_geodesic_upper_bounded_by_cone_path_on_5_taxa():
    from phylomath import cone_path_distance
    t1 = parse_newick('((A:1,B:1):1,(C:1,(D:1,E:1):1):1);')
    t2 = parse_newick('((A:1,C:1):1,(B:1,(D:1,E:1):1):1);')
    d_bhv = bhv_geodesic(t1, t2)
    d_cone = cone_path_distance(t1, t2)
    assert d_bhv <= d_cone + 1e-09

def test_bhv_geodesic_pythagorean_across_orthogonal_orthants():
    t1 = parse_newick('(((A:0,B:0):3,(C:0,D:0):0):1,(E:0,F:0):0);')
    t2 = parse_newick('(((A:0,C:0):4,(B:0,D:0):0):1,(E:0,F:0):0);')
    assert math.isclose(bhv_geodesic(t1, t2), 7.0, rel_tol=1e-09)