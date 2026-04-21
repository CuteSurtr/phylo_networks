from __future__ import annotations
from math import sqrt
from typing import List, Sequence, Tuple
import numpy as np
from .poset import build_closure_lattice, maximal_chains, support_pairs
from .split import Split
from .tree import Tree

def path_space_geo(alpha: Sequence[float], beta: Sequence[float]) -> float:
    return path_space_geo_detail(alpha, beta)[0]

def path_space_geo_detail(alpha: Sequence[float], beta: Sequence[float]) -> Tuple[float, List[List[int]], List[Tuple[float, float]]]:
    a = [float(x) for x in alpha]
    b = [float(x) for x in beta]
    assert len(a) == len(b), 'alpha and beta must have equal length'
    if not a:
        return (0.0, [], [])
    blocks: List[Tuple[float, float, List[int]]] = [(ai, bi, [i]) for i, (ai, bi) in enumerate(zip(a, b))]

    def ratio(pair: Tuple[float, float, List[int]]) -> float:
        A, B, _ = pair
        return float('inf') if B == 0.0 else A / B
    i = 0
    while i + 1 < len(blocks):
        if ratio(blocks[i]) >= ratio(blocks[i + 1]):
            A1, B1, idx1 = blocks[i]
            A2, B2, idx2 = blocks[i + 1]
            merged = (sqrt(A1 * A1 + A2 * A2), sqrt(B1 * B1 + B2 * B2), idx1 + idx2)
            blocks[i:i + 2] = [merged]
            if i > 0:
                i -= 1
        else:
            i += 1
    length = sqrt(sum(((A + B) ** 2 for A, B, _ in blocks)))
    block_indices = [idx for _, _, idx in blocks]
    reduced = [(A, B) for A, B, _ in blocks]
    return (length, block_indices, reduced)

def _split_key(s: Split):
    return frozenset({s.A, s.B})

def geodesic_point_no_common(sigma1: Sequence[Split], sigma2: Sequence[Split], lam: float) -> dict:
    if lam <= 0.0:
        return {_split_key(s): s.length for s in sigma1}
    if lam >= 1.0:
        return {_split_key(s): s.length for s in sigma2}
    if not sigma1 and (not sigma2):
        return {}
    if not sigma1:
        return {_split_key(s): lam * s.length for s in sigma2}
    if not sigma2:
        return {_split_key(s): (1.0 - lam) * s.length for s in sigma1}
    lattice = build_closure_lattice(sigma1, sigma2)
    chains = maximal_chains(lattice)
    best_length = float('inf')
    best_data = None
    for chain in chains:
        pairs = support_pairs(chain, sigma1, sigma2)
        alphas = [_norm_lengths(A) for A, _ in pairs]
        betas = [_norm_lengths(B) for _, B in pairs]
        length, blocks, reduced = path_space_geo_detail(alphas, betas)
        if length < best_length:
            best_length = length
            best_data = (pairs, blocks, reduced)
    if best_data is None:
        if lam < 0.5:
            return {_split_key(s): (1 - 2 * lam) * s.length for s in sigma1}
        return {_split_key(s): (2 * lam - 1) * s.length for s in sigma2}
    pairs, blocks, reduced = best_data
    out: dict = {}
    for j, (A_red, B_red) in enumerate(reduced):
        if A_red + B_red == 0:
            continue
        t_j = A_red / (A_red + B_red)
        for i in blocks[j]:
            A_i, B_i = pairs[i]
            if lam < t_j:
                scale = (t_j - lam) / t_j if t_j > 0 else 0.0
                for e in A_i:
                    out[_split_key(e)] = scale * e.length
            elif lam > t_j:
                scale = (lam - t_j) / (1.0 - t_j) if t_j < 1 else 0.0
                for f in B_i:
                    out[_split_key(f)] = scale * f.length
    return out

def geodesic_point(t1: Tree, t2: Tree, lam: float) -> dict:
    if lam <= 0.0:
        return {_split_key(s): s.length for s in t1.splits()}
    if lam >= 1.0:
        return {_split_key(s): s.length for s in t2.splits()}
    s1 = t1.splits()
    s2 = t2.splits()
    k1 = {_split_key(s): s for s in s1}
    k2 = {_split_key(s): s for s in s2}
    shared_keys = set(k1) & set(k2)
    out: dict = {}
    for k in shared_keys:
        v = (1 - lam) * k1[k].length + lam * k2[k].length
        if v > 0:
            out[k] = v
    for r1, r2 in _residual_blocks(s1, s2):
        for k, v in geodesic_point_no_common(r1, r2, lam).items():
            if v > 0:
                out[k] = v
    return out

def _norm_lengths(splits: Sequence[Split]) -> float:
    return float(np.sqrt(sum((s.length ** 2 for s in splits))))

def geodesic_no_common(sigma1: Sequence[Split], sigma2: Sequence[Split]) -> float:
    if not sigma1 and (not sigma2):
        return 0.0
    if not sigma1:
        return _norm_lengths(sigma2)
    if not sigma2:
        return _norm_lengths(sigma1)
    lattice = build_closure_lattice(sigma1, sigma2)
    chains = maximal_chains(lattice)
    if not chains:
        return _norm_lengths(sigma1) + _norm_lengths(sigma2)
    best = np.inf
    for chain in chains:
        pairs = support_pairs(chain, sigma1, sigma2)
        alpha = [_norm_lengths(A) for A, _ in pairs]
        beta = [_norm_lengths(B) for _, B in pairs]
        L = path_space_geo(alpha, beta)
        best = min(best, L)
    return float(best)

def _residual_blocks(s1: Sequence[Split], s2: Sequence[Split]) -> List[Tuple[List[Split], List[Split]]]:
    k1 = {_split_key(s): s for s in s1}
    k2 = {_split_key(s): s for s in s2}
    shared_keys = set(k1) & set(k2)
    shared = [k1[k] for k in shared_keys]
    residual1 = [s for s in s1 if _split_key(s) not in shared_keys]
    residual2 = [s for s in s2 if _split_key(s) not in shared_keys]

    def small(s: Split):
        return s.A

    def bucket_key(s: Split):
        resid_small = small(s)
        candidates = [sh for sh in shared if resid_small <= small(sh)]
        if not candidates:
            return frozenset()
        best = min(candidates, key=lambda sh: len(small(sh)))
        return small(best)
    buckets1: dict = {}
    buckets2: dict = {}
    for s in residual1:
        buckets1.setdefault(bucket_key(s), []).append(s)
    for s in residual2:
        buckets2.setdefault(bucket_key(s), []).append(s)
    keys = set(buckets1) | set(buckets2)
    return [(buckets1.get(k, []), buckets2.get(k, [])) for k in keys]

def bhv_geodesic(t1: Tree, t2: Tree) -> float:
    s1 = t1.splits()
    s2 = t2.splits()
    k1 = {_split_key(s): s for s in s1}
    k2 = {_split_key(s): s for s in s2}
    shared_keys = set(k1) & set(k2)
    shared_sq = sum(((k1[k].length - k2[k].length) ** 2 for k in shared_keys))
    blocks = _residual_blocks(s1, s2)
    essential_sq = sum((geodesic_no_common(r1, r2) ** 2 for r1, r2 in blocks))
    return float(sqrt(shared_sq + essential_sq))