from __future__ import annotations
from dataclasses import dataclass
from typing import Dict, FrozenSet, Iterable, List, Sequence, Set, Tuple
from .split import Split
SplitKey = FrozenSet[FrozenSet[int]]

def _key(s: Split) -> SplitKey:
    return frozenset({s.A, s.B})

def crossing_set(f: Split, sigma1: Iterable[Split]) -> FrozenSet[SplitKey]:
    return frozenset((_key(e) for e in sigma1 if not f.compatible_with(e)))

def crossing_set_of(fs: Iterable[Split], sigma1: Iterable[Split]) -> FrozenSet[SplitKey]:
    fs = list(fs)
    sigma1 = list(sigma1)
    out: Set[SplitKey] = set()
    for f in fs:
        out |= crossing_set(f, sigma1)
    return frozenset(out)

@dataclass
class ClosureLattice:
    sigma1: List[Split]
    sigma2: List[Split]
    elements: List[FrozenSet[SplitKey]]
    crossing: Dict[FrozenSet[SplitKey], FrozenSet[SplitKey]]
    covers: Dict[int, List[int]]

    @property
    def index_of(self) -> Dict[FrozenSet[SplitKey], int]:
        return {e: i for i, e in enumerate(self.elements)}

def build_closure_lattice(sigma1: Sequence[Split], sigma2: Sequence[Split]) -> ClosureLattice:
    s1 = list(sigma1)
    s2 = list(sigma2)
    keys2 = [_key(f) for f in s2]
    by_key: Dict[SplitKey, Split] = {k: f for k, f in zip(keys2, s2)}

    def closure(A: FrozenSet[SplitKey]) -> FrozenSet[SplitKey]:
        X_A = crossing_set_of((by_key[k] for k in A), s1)
        return frozenset((k for k in keys2 if crossing_set(by_key[k], s1) <= X_A))
    seen: Set[FrozenSet[SplitKey]] = set()
    empty = frozenset()
    frontier = [closure(empty)]
    seen.add(frontier[0])
    while frontier:
        A = frontier.pop()
        for k in keys2:
            if k in A:
                continue
            B = closure(A | {k})
            if B not in seen:
                seen.add(B)
                frontier.append(B)
    elements = sorted(seen, key=lambda s: (len(s), sorted(map(id, s))))
    crossing = {A: crossing_set_of((by_key[k] for k in A), s1) for A in elements}
    covers: Dict[int, List[int]] = {i: [] for i in range(len(elements))}
    for i, A in enumerate(elements):
        for j, B in enumerate(elements):
            if A >= B or not A < B:
                continue
            between = any((A < C < B for C in elements if C is not A and C is not B))
            if not between:
                covers[i].append(j)
    return ClosureLattice(sigma1=s1, sigma2=s2, elements=elements, crossing=crossing, covers=covers)

def maximal_chains(lattice: ClosureLattice) -> List[List[FrozenSet[SplitKey]]]:
    n = len(lattice.elements)
    idx_empty = next((i for i, e in enumerate(lattice.elements) if len(e) == 0))
    full = frozenset((_key(f) for f in lattice.sigma2))
    idx_full = next((i for i, e in enumerate(lattice.elements) if e == full))
    chains: List[List[FrozenSet[SplitKey]]] = []

    def dfs(at: int, acc: List[int]) -> None:
        if at == idx_full:
            chains.append([lattice.elements[i] for i in acc])
            return
        for nxt in lattice.covers[at]:
            dfs(nxt, acc + [nxt])
    dfs(idx_empty, [idx_empty])
    return chains

def support_pairs(chain: Sequence[FrozenSet[SplitKey]], sigma1: Sequence[Split], sigma2: Sequence[Split]) -> List[Tuple[List[Split], List[Split]]]:
    keys2 = [_key(f) for f in sigma2]
    by_key2 = {_key(f): f for f in sigma2}
    by_key1 = {_key(e): e for e in sigma1}
    pairs: List[Tuple[List[Split], List[Split]]] = []
    for i in range(1, len(chain)):
        B_keys = chain[i] - chain[i - 1]
        X_prev = crossing_set_of((by_key2[k] for k in chain[i - 1]), sigma1)
        X_cur = crossing_set_of((by_key2[k] for k in chain[i]), sigma1)
        A_keys = X_cur - X_prev
        A = [by_key1[k] for k in A_keys]
        B = [by_key2[k] for k in B_keys]
        pairs.append((A, B))
    return pairs