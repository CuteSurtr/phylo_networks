from __future__ import annotations
import random
from typing import Dict, FrozenSet, List, Sequence
from .gtp import geodesic_point
from .split import Split
from .tree import Tree
SplitKey = FrozenSet[FrozenSet[int]]

def _key(s: Split) -> SplitKey:
    return frozenset({s.A, s.B})

def _lengths_to_virtual_tree(lengths: Dict[SplitKey, float]) -> '_LengthMapTree':
    return _LengthMapTree(lengths)

class _LengthMapTree:

    def __init__(self, lengths: Dict[SplitKey, float]):
        self._lengths = {k: float(v) for k, v in lengths.items() if v > 0}

    def splits(self) -> List[Split]:
        out: List[Split] = []
        for k, v in self._lengths.items():
            A, B = tuple(k)
            out.append(Split(A=A, B=B, length=v))
        return out

def frechet_mean(trees: Sequence[Tree], n_iters: int=500, seed: int | None=0) -> Dict[SplitKey, float]:
    if not trees:
        raise ValueError('need at least one tree')
    N = len(trees)
    rng = random.Random(seed)
    start = trees[rng.randrange(N)]
    mu_lengths: Dict[SplitKey, float] = {_key(s): s.length for s in start.splits()}
    for k in range(n_iters):
        target = trees[(k + 1 + rng.randrange(max(1, N))) % N]
        step = 1.0 / (k + 2)
        mu_tree = _lengths_to_virtual_tree(mu_lengths)
        new_lengths = geodesic_point(mu_tree, target, step)
        mu_lengths = {k: v for k, v in new_lengths.items() if v > 1e-10}
    return mu_lengths

def frechet_variance(mean_lengths: Dict[SplitKey, float], trees: Sequence[Tree]) -> float:
    if not trees:
        raise ValueError('need at least one tree')
    from .gtp import bhv_geodesic
    mu_tree = _lengths_to_virtual_tree(mean_lengths)
    total = 0.0
    for t in trees:
        d = bhv_geodesic(mu_tree, t)
        total += d * d
    return total / len(trees)