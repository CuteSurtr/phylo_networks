from __future__ import annotations
from typing import Dict, List, Optional, Sequence, Tuple
import matplotlib.pyplot as plt
import numpy as np
from .tree import Node, Tree

def _leaf_x_coords(t: Tree) -> Dict[int, float]:
    order: List[Node] = []

    def walk(n: Node) -> None:
        if n.is_leaf:
            order.append(n)
        else:
            for c in n.children:
                walk(c)
    walk(t.root)
    return {id(lf): float(i) for i, lf in enumerate(order)}

def _node_depths(t: Tree) -> Dict[int, float]:
    depth: Dict[int, float] = {id(t.root): 0.0}
    for n in t.postorder():
        pass
    stack = [t.root]
    while stack:
        n = stack.pop()
        for c in n.children:
            depth[id(c)] = depth[id(n)] + max(c.length, 1e-06)
            stack.append(c)
    return depth

def plot_tree(t: Tree, ax: Optional[plt.Axes]=None, label_leaves: bool=True):
    if ax is None:
        _, ax = plt.subplots(figsize=(4.5, 4))
    leaf_x = _leaf_x_coords(t)
    depths = _node_depths(t)
    max_d = max(depths.values()) if depths else 1.0
    node_x: Dict[int, float] = {}

    def assign(n: Node) -> float:
        if n.is_leaf:
            node_x[id(n)] = leaf_x[id(n)]
            return node_x[id(n)]
        xs = [assign(c) for c in n.children]
        node_x[id(n)] = sum(xs) / len(xs)
        return node_x[id(n)]
    assign(t.root)

    def draw(n: Node) -> None:
        x_p = node_x[id(n)]
        y_p = -depths[id(n)]
        for c in n.children:
            x_c = node_x[id(c)]
            y_c = -depths[id(c)]
            ax.plot([x_c, x_c], [y_p, y_c], color='black', linewidth=1.0)
            ax.plot([x_p, x_c], [y_p, y_p], color='black', linewidth=1.0)
            draw(c)
    draw(t.root)
    if label_leaves:
        for n in t.postorder():
            if n.is_leaf:
                ax.text(node_x[id(n)], -max_d - 0.05 * max_d - 0.1, n.name, ha='center', va='top', fontsize=9)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_ylim(-max_d * 1.25, 0.1)
    return ax

def plot_tanglegram(t1: Tree, t2: Tree, ax: Optional[plt.Axes]=None) -> plt.Axes:
    from .tanglegram import leaf_order
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 5))
    depths1 = _node_depths(t1)
    depths2 = _node_depths(t2)
    max_d1 = max(depths1.values(), default=1.0)
    max_d2 = max(depths2.values(), default=1.0)
    order1 = leaf_order(t1)
    order2 = leaf_order(t2)
    n = max(len(order1), len(order2))
    y_of_leaf_1 = {tid: i for i, tid in enumerate(order1)}
    y_of_leaf_2 = {tid: i for i, tid in enumerate(order2)}

    def render(t: Tree, x_leaf: float, x_root: float, y_of_leaf: Dict[int, float]) -> Dict[int, Tuple[float, float]]:
        depths = _node_depths(t)
        max_d = max(depths.values(), default=1.0)
        coord: Dict[int, Tuple[float, float]] = {}

        def assign(n: Node) -> Tuple[float, float]:
            if n.is_leaf:
                xy = (x_leaf, y_of_leaf[n.taxon_id])
                coord[id(n)] = xy
                return xy
            xys = [assign(c) for c in n.children]
            y = sum((xy[1] for xy in xys)) / len(xys)
            frac = depths[id(n)] / max_d
            x = x_leaf + (x_root - x_leaf) * (1.0 - frac)
            coord[id(n)] = (x, y)
            return coord[id(n)]
        assign(t.root)

        def draw(n: Node):
            x_p, y_p = coord[id(n)]
            for c in n.children:
                x_c, y_c = coord[id(c)]
                ax.plot([x_p, x_p], [y_p, y_c], color='black', linewidth=1.0)
                ax.plot([x_p, x_c], [y_c, y_c], color='black', linewidth=1.0)
                draw(c)
        draw(t.root)
        return coord
    render(t1, x_leaf=0.35, x_root=0.0, y_of_leaf=y_of_leaf_1)
    render(t2, x_leaf=0.65, x_root=1.0, y_of_leaf=y_of_leaf_2)
    for tid, y in y_of_leaf_1.items():
        ax.text(0.33, y, t1.taxa[tid], ha='right', va='center', fontsize=9)
    for tid, y in y_of_leaf_2.items():
        ax.text(0.67, y, t2.taxa[tid], ha='left', va='center', fontsize=9)
    crossings = 0
    for tid in set(y_of_leaf_1) & set(y_of_leaf_2):
        y1 = y_of_leaf_1[tid]
        y2 = y_of_leaf_2[tid]
        ax.plot([0.35, 0.65], [y1, y2], color='tab:blue', alpha=0.5, linewidth=1.0)
    from .tanglegram import count_crossings
    crossings = count_crossings(order1, order2)
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.8, n - 0.2)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(f'Tanglegram (crossings: {crossings})')
    return ax

def plot_bhv_T4_link(ax: Optional[plt.Axes]=None) -> plt.Axes:
    if ax is None:
        _, ax = plt.subplots(figsize=(5, 5))
    outer = [(np.cos(a + np.pi / 2), np.sin(a + np.pi / 2)) for a in np.linspace(0, 2 * np.pi, 6)[:-1]]
    inner = [(0.55 * np.cos(a + np.pi / 2), 0.55 * np.sin(a + np.pi / 2)) for a in np.linspace(0, 2 * np.pi, 6)[:-1]]
    for i in range(5):
        x1, y1 = outer[i]
        x2, y2 = outer[(i + 1) % 5]
        ax.plot([x1, x2], [y1, y2], color='black', linewidth=1.3)
    for i in range(5):
        x1, y1 = inner[i]
        x2, y2 = inner[(i + 2) % 5]
        ax.plot([x1, x2], [y1, y2], color='black', linewidth=1.3)
    for i in range(5):
        ax.plot([outer[i][0], inner[i][0]], [outer[i][1], inner[i][1]], color='black', linewidth=1.3)
    labels_outer = ['{1,2}', '{2,3}', '{3,4}', '{4,5}', '{1,5}']
    labels_inner = ['{3,4,5}', '{1,4,5}', '{1,2,5}', '{1,2,3}', '{2,3,4}']
    for (x, y), lab in zip(outer, labels_outer):
        ax.plot(x, y, 'o', color='tab:red', markersize=8)
        ax.text(x * 1.15, y * 1.15, lab, ha='center', va='center', fontsize=9)
    for (x, y), lab in zip(inner, labels_inner):
        ax.plot(x, y, 'o', color='tab:blue', markersize=8)
        ax.text(x * 0.75, y * 0.75, lab, ha='center', va='center', fontsize=8)
    ax.set_xlim(-1.6, 1.6)
    ax.set_ylim(-1.4, 1.4)
    ax.set_aspect('equal')
    ax.set_title('Link of origin in T? = Petersen graph')
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    return ax

def plot_frechet_path(trees: Sequence[Tree], n_iters: int=200, seed: int=0, ax: Optional[plt.Axes]=None) -> plt.Axes:
    import random
    from .frechet import _LengthMapTree, _key
    from .gtp import bhv_geodesic, geodesic_point
    if ax is None:
        _, ax = plt.subplots(figsize=(5, 3.5))
    if not trees:
        raise ValueError('need at least one tree')
    N = len(trees)
    rng = random.Random(seed)
    start = trees[rng.randrange(N)]
    mu_lengths: Dict = {_key(s): s.length for s in start.splits()}

    def obj(mu_lengths) -> float:
        mu = _LengthMapTree(mu_lengths)
        return sum((bhv_geodesic(mu, t) ** 2 for t in trees)) / N
    traj = [obj(mu_lengths)]
    for k in range(n_iters):
        target = trees[(k + 1 + rng.randrange(N)) % N]
        step = 1.0 / (k + 2)
        mu_tree = _LengthMapTree(mu_lengths)
        mu_lengths = geodesic_point(mu_tree, target, step)
        traj.append(obj(mu_lengths))
    ax.plot(traj, color='tab:purple', linewidth=1.3)
    ax.set_xlabel('iteration k')
    ax.set_ylabel('Frechet objective')
    ax.set_title('Bacak proximal-point convergence')
    ax.grid(alpha=0.3)
    return ax

def plot_nni_graph(n: int, ax: Optional[plt.Axes]=None) -> plt.Axes:
    try:
        import networkx as nx
    except ImportError as e:
        raise RuntimeError('networkx required for plot_nni_graph') from e
    from .rearrange import nni_graph
    G = nni_graph(n)
    if ax is None:
        _, ax = plt.subplots(figsize=(6, 6))
    pos = nx.spring_layout(G, seed=0)
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color='gray', alpha=0.4, width=0.6)
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=20, node_color='tab:red')
    ax.set_title(f'NNI graph, n={n}: {G.number_of_nodes()} vertices, {G.number_of_edges()} edges')
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    return ax

def save(fig_or_ax, path: str) -> None:
    if hasattr(fig_or_ax, 'savefig'):
        fig_or_ax.savefig(path, dpi=150, bbox_inches='tight')
    else:
        fig_or_ax.figure.savefig(path, dpi=150, bbox_inches='tight')