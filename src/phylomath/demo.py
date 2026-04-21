from __future__ import annotations
import os
from math import sqrt
from pathlib import Path
import numpy as np
from . import bhv_distance, bhv_geodesic, binary_topologies, cone_path_distance, count_crossings, fermat_weber, frechet_mean, frechet_variance, leaf_order, min_crossings_exhaustive, min_crossings_one_tree, neighbor_net, nni_graph, parse_newick, path_space_geo, quartet_distance, rf_distance, to_newick, tree_dissimilarity, tropical_distance, weighted_rf_distance
from .split import pairwise_compatible
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from . import viz
    MATPLOTLIB = True
except ImportError:
    MATPLOTLIB = False
HERE = Path(__file__).resolve().parents[2]
DATA = HERE / 'data'
RESULTS = HERE / 'results'

def section(title: str) -> None:
    print('\n' + '=' * 74)
    print(title)
    print('=' * 74)

def main() -> None:
    RESULTS.mkdir(exist_ok=True)
    section('1. Newick + splits + compatibility (flag-complex property)')
    t1 = parse_newick('((A:1.0,B:1.0):0.5,(C:2.0,(D:0.1,E:0.3):0.2):0.4);')
    t2 = parse_newick('((A:1.1,C:0.9):0.6,(B:2.1,(D:0.2,E:0.4):0.3):0.5);')
    print('t1 =', to_newick(t1))
    print('t2 =', to_newick(t2))
    print('pairwise compatible?', pairwise_compatible(t1.splits()))
    section('2. Tree distances (all metrics)')
    print(f'RF              = {rf_distance(t1, t2)}')
    print(f'wRF  (l2)       = {weighted_rf_distance(t1, t2):.4f}')
    print(f'quartet         = {quartet_distance(t1, t2)}')
    print(f'cone path       = {cone_path_distance(t1, t2):.4f}')
    print(f'BHV (UB, old)   = {bhv_distance(t1, t2):.4f}')
    print(f'BHV (GTP exact) = {bhv_geodesic(t1, t2):.4f}')
    section('3. PathSpaceGeo closed-form sanity')
    print(f'k=1 through origin          : α=3,β=4 → {path_space_geo([3.0], [4.0])}')
    print(f'k=2 ascending ratios        : (1,3)/(2,1) → {path_space_geo([1.0, 3.0], [2.0, 1.0])}')
    print(f'k=2 descending (merges)     : (3,1)/(1,2) → {path_space_geo([3.0, 1.0], [1.0, 2.0]):.4f}')
    section('4. Enumerate rooted binary trees (Schroder counting)')
    for n in (3, 4, 5, 6, 7):
        c = sum((1 for _ in binary_topologies(n)))
        print(f'  n={n}: {c}')
    section('5. NNI rearrangement graph on n=5')
    G = nni_graph(5)
    print(f'  {G.number_of_nodes()} vertices, {G.number_of_edges()} edges')
    section('6. Frechet mean + variance on primate mitochondrial sample')
    trees = []
    for line in (DATA / 'primates_mitochondrial.nwk').read_text().splitlines():
        if line.strip():
            trees.append(parse_newick(line.strip()))
    print(f'  loaded {len(trees)} trees from data/primates_mitochondrial.nwk')
    mu = frechet_mean(trees, n_iters=800, seed=0)
    var = frechet_variance(mu, trees)
    print(f'  Frechet variance = {var:.6f}')
    print(f'  mean has {len(mu)} non-zero internal splits')
    section('7. Tanglegram — default vs OTCM vs TTCM')
    h = parse_newick('(((A,B),C),(D,E));')
    p = parse_newick('(((A,C),B),(E,D));')
    print(f'  default layout crossings   = {count_crossings(leaf_order(h), leaf_order(p))}')
    print(f'  OTCM DP (free=parasite)    = {min_crossings_one_tree(h, p)}')
    print(f'  TTCM exhaustive            = {min_crossings_exhaustive(h, p)}')
    section('8. Neighbor-Net circular ordering')
    u = tree_dissimilarity(trees[0])
    n = len(trees[0].taxa)
    D = np.zeros((n, n))
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            D[i, j] = D[j, i] = u[k]
            k += 1
    order, splits = neighbor_net(D)
    labels = [trees[0].taxa[i] for i in order]
    print(f'  circular order   : {labels}')
    print(f'  circular splits  : {len(splits)}')
    section('9. Tropical geometry + Fermat-Weber')
    Us = [tree_dissimilarity(t) for t in trees]
    d_trop = tropical_distance(Us[0], Us[1])
    fw = fermat_weber(Us)
    obj_fw = sum((tropical_distance(fw, u) for u in Us))
    obj_u0 = sum((tropical_distance(Us[0], u) for u in Us))
    print(f'  d∞(T_0, T_1)        = {d_trop:.4f}')
    print(f'  FW objective        = {obj_fw:.4f}')
    print(f'  objective at T_0    = {obj_u0:.4f}   (FW should be ≤)')
    assert obj_fw <= obj_u0 + 1e-06, "FW point isn't reducing objective"
    if MATPLOTLIB:
        section('10. Visualization — writing figures to results/')
        fig, axs = plt.subplots(2, 2, figsize=(11, 9))
        viz.plot_tree(trees[0], ax=axs[0][0])
        axs[0][0].set_title('Primate tree 1')
        viz.plot_tanglegram(h, p, ax=axs[0][1])
        viz.plot_bhv_T4_link(ax=axs[1][0])
        viz.plot_frechet_path(trees, n_iters=150, ax=axs[1][1])
        fig.tight_layout()
        fig.savefig(RESULTS / 'phylomath_demo.png', dpi=150, bbox_inches='tight')
        print(f"  wrote {RESULTS / 'phylomath_demo.png'}")
        fig2, ax2 = plt.subplots(figsize=(6, 6))
        viz.plot_nni_graph(5, ax=ax2)
        fig2.savefig(RESULTS / 'nni_graph_n5.png', dpi=150, bbox_inches='tight')
        print(f"  wrote {RESULTS / 'nni_graph_n5.png'}")
    else:
        print('\n  (matplotlib not available — skipping viz)')
if __name__ == '__main__':
    main()