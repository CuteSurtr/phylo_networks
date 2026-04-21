from __future__ import annotations
from typing import Dict, List, Sequence, Tuple
import numpy as np
from scipy.optimize import linprog
from .tree import Node, Tree

def tree_dissimilarity(t: Tree) -> np.ndarray:
    n = len(t.taxa)
    adj: Dict[int, List[Tuple[int, float]]] = {id(n_): [] for n_ in t.postorder()}
    for n_ in t.postorder():
        for c in n_.children:
            adj[id(n_)].append((id(c), c.length))
            adj[id(c)].append((id(n_), c.length))
    leaf_nodes = t.leaves()
    leaf_by_tid = {lf.taxon_id: id(lf) for lf in leaf_nodes}
    m = n * (n - 1) // 2
    u = np.zeros(m)

    def idx(i: int, j: int) -> int:
        a, b = (i, j) if i < j else (j, i)
        return a * n - a * (a + 1) // 2 + (b - a - 1)
    for i in range(n):
        start = leaf_by_tid[i]
        dist = {start: 0.0}
        stack = [start]
        while stack:
            x = stack.pop()
            for y, w in adj[x]:
                if y in dist:
                    continue
                dist[y] = dist[x] + w
                stack.append(y)
        for j in range(i + 1, n):
            u[idx(i, j)] = dist[leaf_by_tid[j]]
    return u

def tropical_distance(u: np.ndarray, v: np.ndarray) -> float:
    u = np.asarray(u, dtype=float)
    v = np.asarray(v, dtype=float)
    assert u.shape == v.shape
    diff = u - v
    return float(diff.max() - diff.min())

def fermat_weber(points: Sequence[np.ndarray]) -> np.ndarray:
    P = np.stack([np.asarray(p, dtype=float) for p in points], axis=0)
    N, m = P.shape
    nvar = m + 2 * N
    c = np.zeros(nvar)
    for i in range(N):
        c[m + i] = 1.0
        c[m + N + i] = -1.0
    rows_A = []
    rows_b = []
    for i in range(N):
        for j in range(m):
            row = np.zeros(nvar)
            row[j] = 1.0
            row[m + i] = -1.0
            rows_A.append(row)
            rows_b.append(P[i, j])
            row = np.zeros(nvar)
            row[j] = -1.0
            row[m + N + i] = 1.0
            rows_A.append(row)
            rows_b.append(-P[i, j])
    A_ub = np.asarray(rows_A)
    b_ub = np.asarray(rows_b)
    A_eq = np.zeros((1, nvar))
    A_eq[0, 0] = 1.0
    b_eq = np.array([0.0])
    bounds = [(None, None)] * nvar
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method='highs')
    if not res.success:
        raise RuntimeError(f'FW LP failed: {res.message}')
    return res.x[:m]

def tropical_line_through(u: np.ndarray, v: np.ndarray, n_samples: int=50) -> np.ndarray:
    u = np.asarray(u, dtype=float)
    v = np.asarray(v, dtype=float)
    d = tropical_distance(u, v)
    params = np.linspace(-d, d, n_samples)
    samples = np.stack([np.minimum(u + t, v) for t in params], axis=0)
    return samples