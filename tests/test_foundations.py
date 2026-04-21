import math
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from phylomath import binary_topologies, bhv_distance, cone_path_distance, nni_neighbors, parse_newick, quartet_distance, rf_distance, to_newick, weighted_rf_distance
from phylomath.split import Split, pairwise_compatible

def test_newick_roundtrip():
    s = '((A:1,B:2):0.5,(C:1.5,D:0.7):0.3);'
    t = parse_newick(s)
    assert len(t.leaves()) == 4
    t2 = parse_newick(to_newick(t))
    assert sorted((n.name for n in t.leaves())) == sorted((n.name for n in t2.leaves()))

def test_splits_count():
    t = parse_newick('((A,B),(C,(D,E)));')
    assert len(t.taxa) == 5
    assert len([s for s in t.splits() if not s.is_trivial]) == 5 - 3

def test_compatibility():
    taxa = list(range(5))
    s1 = Split.of({0, 1}, taxa)
    s2 = Split.of({0, 1, 2}, taxa)
    assert s1.compatible_with(s2)
    s3 = Split.of({0, 2}, taxa)
    assert not s1.compatible_with(s3)

def test_binary_topology_count():
    assert sum((1 for _ in binary_topologies(3))) == 3
    assert sum((1 for _ in binary_topologies(4))) == 15
    assert sum((1 for _ in binary_topologies(5))) == 105
    assert sum((1 for _ in binary_topologies(6))) == 945

def test_pairwise_compatibility_of_generated_trees():
    for s in binary_topologies(5):
        t = parse_newick(s)
        assert pairwise_compatible(t.splits())

def test_rf_identity_and_symmetry():
    t1 = parse_newick('((A,B),(C,(D,E)));')
    t2 = parse_newick('((A,B),(C,(D,E)));')
    assert rf_distance(t1, t2) == 0
    t3 = parse_newick('((A,C),(B,(D,E)));')
    assert rf_distance(t1, t3) == rf_distance(t3, t1) > 0

def test_weighted_rf_zero_on_identity():
    t = parse_newick('((A:1,B:2):0.5,(C:3,D:0.7):0.2);')
    assert weighted_rf_distance(t, t) == 0.0
    assert weighted_rf_distance(t, t, 'l1') == 0.0

def test_quartet_distance_zero_on_identity():
    t = parse_newick('((A,B),(C,(D,E)));')
    assert quartet_distance(t, t) == 0

def test_cone_path_zero_on_identical_tree():
    t = parse_newick('((A:1,B:2):0.5,(C:3,D:0.7):0.2);')
    assert math.isclose(cone_path_distance(t, t), 0.0, abs_tol=1e-09)

def test_bhv_upper_bounded_by_cone():
    t1 = parse_newick('((A:1,B:1):0.5,(C:2,(D:0.1,E:0.3):0.2):0.4);')
    t2 = parse_newick('((A:1.1,C:0.9):0.6,(B:2.1,(D:0.2,E:0.4):0.3):0.5);')
    d_cone = cone_path_distance(t1, t2)
    d_bhv = bhv_distance(t1, t2)
    assert d_bhv <= d_cone + 1e-09

def test_nni_neighbors_sane_count():
    t = parse_newick('((1,2),(3,(4,5)));')
    nbrs = nni_neighbors(t)
    assert 1 <= len(nbrs) <= 2 * (len(t.taxa) - 2)