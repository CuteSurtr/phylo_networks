import math
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from phylomath import bhv_geodesic, count_crossings, frechet_mean, frechet_variance, geodesic_no_common, leaf_order, min_crossings_exhaustive, min_crossings_one_tree, nni_graph, parse_newick, shortest_nni_distance, to_newick

def test_gtp_zero_on_identical_trees():
    t = parse_newick('((A:1,B:2):0.5,(C:3,D:0.7):0.2);')
    d = bhv_geodesic(t, t)
    assert d == 0.0

def test_gtp_matches_pure_branch_length_on_same_topology():
    t1 = parse_newick('((A:1,B:1):1,(C:1,D:1):1);')
    t2 = parse_newick('((A:2,B:2):2,(C:2,D:2):2);')
    d = bhv_geodesic(t1, t2)
    assert math.isclose(d, 2.0, rel_tol=1e-06)

def test_gtp_cone_path_for_totally_incompatible_tiny_case():
    t1 = parse_newick('((A:0,B:0):3,(C:0,D:0):0);')
    t2 = parse_newick('((A:0,C:0):4,(B:0,D:0):0);')
    d = bhv_geodesic(t1, t2)
    assert math.isclose(d, 7.0, rel_tol=1e-06)

def test_frechet_mean_on_identical_samples():
    t = parse_newick('((A:1,B:2):0.5,(C:3,D:0.7):0.2);')
    mu = frechet_mean([t, t, t], n_iters=200)
    var = frechet_variance(mu, [t, t, t])
    assert math.isclose(var, 0.0, abs_tol=1e-06)

def test_frechet_mean_averages_branch_lengths():
    t1 = parse_newick('((A:1,B:1):2,(C:1,D:1):1);')
    t2 = parse_newick('((A:1,B:1):4,(C:1,D:1):3);')
    mu = frechet_mean([t1, t2], n_iters=2000, seed=42)
    lens = [v for v in mu.values() if v > 0.01]
    assert any((abs(v - 5.0) < 0.3 for v in lens)), f'no length near 5.0 in {lens}'

def test_count_crossings_same_order_zero():
    o = [0, 1, 2, 3, 4]
    assert count_crossings(o, o) == 0

def test_count_crossings_full_reversal():
    n = 5
    assert count_crossings(list(range(n)), list(reversed(range(n)))) == n * (n - 1) // 2

def test_tanglegram_identical_trees_zero_minimum():
    t1 = parse_newick('((A,B),(C,D));')
    t2 = parse_newick('((A,B),(C,D));')
    assert min_crossings_exhaustive(t1, t2) == 0
    assert min_crossings_one_tree(t1, t2) == 0

def test_tanglegram_different_trees_has_positive_crossings():
    t1 = parse_newick('((A,B),(C,D));')
    t2 = parse_newick('((A,C),(B,D));')
    assert min_crossings_exhaustive(t1, t2) >= 1

def test_nni_graph_small_sizes():
    G = nni_graph(3)
    assert G.number_of_nodes() == 3
    G4 = nni_graph(4)
    assert G4.number_of_nodes() == 15
    G5 = nni_graph(5)
    assert G5.number_of_nodes() == 105
    assert G5.number_of_edges() > 0

def test_shortest_nni_distance_zero_self():
    G = nni_graph(4)
    s = '((1,2),(3,4));'
    assert shortest_nni_distance(s, s, G) == 0