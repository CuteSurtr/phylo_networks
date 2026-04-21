from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, FrozenSet, List, Optional, Sequence, Set, Tuple
import numpy as np

@dataclass
class NetNode:
    name: Optional[str] = None
    children: List['NetNode'] = field(default_factory=list)
    parents: List['NetNode'] = field(default_factory=list)

    @property
    def is_leaf(self) -> bool:
        return not self.children

    @property
    def is_reticulation(self) -> bool:
        return len(self.parents) >= 2

@dataclass
class PhyloNetwork:
    root: NetNode
    taxa: List[str]

    def nodes(self) -> List[NetNode]:
        seen: Set[int] = set()
        order: List[NetNode] = []

        def walk(n: NetNode) -> None:
            if id(n) in seen:
                return
            seen.add(id(n))
            order.append(n)
            for c in n.children:
                walk(c)
        walk(self.root)
        return order

    def reticulations(self) -> List[NetNode]:
        return [n for n in self.nodes() if n.is_reticulation]

    def leaves(self) -> List[NetNode]:
        return [n for n in self.nodes() if n.is_leaf]

def reticulation_count(net: PhyloNetwork) -> int:
    return len(net.reticulations())

def is_tree_based(net: PhyloNetwork) -> bool:
    rets = net.reticulations()
    if not rets:
        return True

    def backtrack(i: int, kept_edges: List[Tuple[NetNode, NetNode]]) -> bool:
        if i == len(rets):
            return _is_spanning_tree(net, kept_edges)
        r = rets[i]
        for p in r.parents:
            if backtrack(i + 1, kept_edges + [(p, r)]):
                return True
        return False
    return backtrack(0, [])

def _is_spanning_tree(net: PhyloNetwork, reticulation_choices: List[Tuple[NetNode, NetNode]]) -> bool:
    nodes = net.nodes()
    adj: Dict[int, List[int]] = {id(n): [] for n in nodes}
    chosen_ret = {id(r): id(p) for p, r in reticulation_choices}
    for n in nodes:
        for c in n.children:
            if c.is_reticulation:
                if chosen_ret.get(id(c)) != id(n):
                    continue
            adj[id(n)].append(id(c))
            adj[id(c)].append(id(n))
    start = id(net.root)
    visited: Set[int] = set()
    stack = [start]
    while stack:
        x = stack.pop()
        if x in visited:
            continue
        visited.add(x)
        stack.extend(adj[x])
    if len(visited) != len(nodes):
        return False
    edge_count = sum((len(v) for v in adj.values())) // 2
    return edge_count == len(nodes) - 1

def neighbor_net(D: np.ndarray, labels: Optional[Sequence[str]]=None) -> Tuple[List[int], List[FrozenSet[int]]]:
    D = np.asarray(D, dtype=float).copy()
    n = D.shape[0]
    assert D.shape == (n, n)
    clusters: List[List[int]] = [[i] for i in range(n)]

    def Q_matrix(Dm: np.ndarray) -> np.ndarray:
        m = Dm.shape[0]
        row_sum = Dm.sum(axis=1)
        Q = (m - 2) * Dm - row_sum[:, None] - row_sum[None, :]
        np.fill_diagonal(Q, np.inf)
        return Q
    while D.shape[0] > 3:
        Q = Q_matrix(D)
        i, j = np.unravel_index(np.argmin(Q), Q.shape)
        if i > j:
            i, j = (j, i)
        merged = clusters[i] + clusters[j]
        m = D.shape[0]
        new_dist = np.zeros(m - 1)
        keep = [k for k in range(m) if k != i and k != j]
        for idx, k in enumerate(keep):
            new_dist[idx] = 0.5 * (D[i, k] + D[j, k] - D[i, j])
        D_new = np.zeros((m - 1, m - 1))
        for a, ka in enumerate(keep):
            for b, kb in enumerate(keep):
                D_new[a, b] = D[ka, kb]
            D_new[a, m - 2] = new_dist[a]
            D_new[m - 2, a] = new_dist[a]
        D = D_new
        new_clusters = [clusters[k] for k in keep]
        new_clusters.append(merged)
        clusters = new_clusters
    order: List[int] = []
    for c in clusters:
        order.extend(c)
    splits: List[FrozenSet[int]] = []
    n_full = len(order)
    seen: Set[FrozenSet[int]] = set()
    for length in range(1, n_full):
        for start in range(n_full):
            block = frozenset((order[(start + k) % n_full] for k in range(length)))
            if len(block) > n_full - len(block):
                block = frozenset(range(n_full)) - block
            if 0 < len(block) < n_full and block not in seen:
                seen.add(block)
                splits.append(block)
    return (order, splits)