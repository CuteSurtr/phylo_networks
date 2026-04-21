import math
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
import numpy as np
import pytest
from phylomath import PhyloNetwork, NetNode, fermat_weber, is_tree_based, neighbor_net, parse_newick, reticulation_count, tree_dissimilarity, tropical_distance, tropical_line_through
from phylomath.gtp import geodesic_point

def _toy_network_with_one_reticulation() -> PhyloNetwork:
    a = NetNode(name='A')
    b = NetNode(name='B')
    c = NetNode(name='C')
    internal = NetNode(name=None)
    recomb = NetNode(name=None)
    internal.children = [b, recomb]
    b.parents = [internal]
    recomb.parents = [internal]
    recomb.children = [c]
    c.parents = [recomb]
    root = NetNode(name=None)
    root.children = [a, internal, recomb]
    a.parents = [root]
    internal.parents = [root]
    recomb.parents.append(root)
    net = PhyloNetwork(root=root, taxa=['A', 'B', 'C'])
    return net

def test_reticulation_count():
    net = _toy_network_with_one_reticulation()
    assert reticulation_count(net) == 1

def test_is_tree_based_trivial_tree():
    a = NetNode(name='A')
    b = NetNode(name='B')
    root = NetNode(children=[a, b])
    a.parents = [root]
    b.parents = [root]
    net = PhyloNetwork(root=root, taxa=['A', 'B'])
    assert is_tree_based(net)
    assert reticulation_count(net) == 0

def test_neighbor_net_circular_matrix_recovers_order():
    n = 6
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                d = min((i - j) % n, (j - i) % n)
                D[i, j] = float(d)
    order, splits = neighbor_net(D)
    assert sorted(order) == list(range(n))
    assert len(splits) > 0

def test_tree_dissimilarity_triangle_inequality():
    t = parse_newick('((A:1,B:2):0.5,(C:1.5,D:0.7):0.3);')
    u = tree_dissimilarity(t)
    assert (u >= 0).all()
    n = len(t.taxa)

    def pair_idx(i, j):
        a, b = (i, j) if i < j else (j, i)
        return a * n - a * (a + 1) // 2 + (b - a - 1)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            D[i, j] = D[j, i] = u[pair_idx(i, j)]
    a, b, c, d = (0, 1, 2, 3)
    s1 = D[a, b] + D[c, d]
    s2 = D[a, c] + D[b, d]
    s3 = D[a, d] + D[b, c]
    sums = sorted([s1, s2, s3])
    assert math.isclose(sums[1], sums[2], rel_tol=1e-09)

def test_tropical_distance_zero_on_identity():
    u = np.array([1.0, 2.0, 3.0, 1.5, 2.5, 0.8])
    assert tropical_distance(u, u) == 0.0

def test_tropical_distance_invariant_to_adding_constant():
    u = np.array([1.0, 2.0, 3.0])
    v = np.array([5.0, 6.0, 7.0])
    assert tropical_distance(u, v) == 0.0

def test_tropical_distance_symmetric():
    u = np.array([0.0, 1.0, 2.0, 3.0])
    v = np.array([1.0, 0.0, 3.0, 1.0])
    assert tropical_distance(u, v) == tropical_distance(v, u)

def test_fermat_weber_recovers_centroid_on_colinear_points():
    u = np.array([0.0, 0.0, 0.0])
    v = np.array([1.0, 2.0, 3.0])
    w = np.array([2.0, 4.0, 6.0])
    fw = fermat_weber([u, v, w])
    obj_fw = sum((tropical_distance(fw, p) for p in (u, v, w)))
    obj_u = sum((tropical_distance(u, p) for p in (u, v, w)))
    obj_v = sum((tropical_distance(v, p) for p in (u, v, w)))
    obj_w = sum((tropical_distance(w, p) for p in (u, v, w)))
    assert obj_fw <= obj_u + 1e-09
    assert obj_fw <= obj_v + 1e-09
    assert obj_fw <= obj_w + 1e-09

def test_tropical_line_sample_shape():
    u = np.array([0.0, 1.0, 2.0])
    v = np.array([1.0, 0.0, 3.0])
    samples = tropical_line_through(u, v, n_samples=20)
    assert samples.shape == (20, 3)

def test_geodesic_point_endpoints():
    t1 = parse_newick('((A:1,B:1):2,(C:1,D:1):1);')
    t2 = parse_newick('((A:1,B:1):4,(C:1,D:1):3);')
    at0 = geodesic_point(t1, t2, 0.0)
    at1 = geodesic_point(t1, t2, 1.0)
    from phylomath.gtp import _split_key
    exp0 = {_split_key(s): s.length for s in t1.splits()}
    exp1 = {_split_key(s): s.length for s in t2.splits()}
    assert at0 == exp0
    assert at1 == exp1

def test_geodesic_point_midpoint_shared_topology():
    t1 = parse_newick('((A:1,B:1):2,(C:1,D:1):1);')
    t2 = parse_newick('((A:1,B:1):4,(C:1,D:1):3);')
    mid = geodesic_point(t1, t2, 0.5)
    vals = list(mid.values())
    assert any((abs(v - 5.0) < 1e-06 for v in vals))

def test_parse_newick_empty_raises():
    with pytest.raises(ValueError):
        parse_newick('')

def test_parse_newick_non_string_raises():
    with pytest.raises(TypeError):
        parse_newick(42)

def test_rf_distance_mismatched_taxa_raises():
    from phylomath import rf_distance
    t1 = parse_newick('((A,B),(C,D));')
    t2 = parse_newick('((A,B),(C,E));')
    with pytest.raises(ValueError):
        rf_distance(t1, t2)