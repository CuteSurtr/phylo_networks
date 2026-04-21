from __future__ import annotations
from itertools import combinations
from math import sqrt
from typing import Literal
from .split import Split
from .tree import Tree

def rf_distance(t1: Tree, t2: Tree) -> int:
    if sorted(t1.taxa) != sorted(t2.taxa):
        raise ValueError(f'rf_distance requires identical taxon sets; got {t1.taxa} vs {t2.taxa}')
    s1 = {_key(s) for s in t1.splits()}
    s2 = {_key(s) for s in t2.splits()}
    return len(s1 ^ s2)

def weighted_rf_distance(t1: Tree, t2: Tree, norm: Literal['l1', 'l2']='l2') -> float:
    m1 = {_key(s): s.length for s in t1.splits()}
    m2 = {_key(s): s.length for s in t2.splits()}
    keys = set(m1) | set(m2)
    diffs = [m1.get(k, 0.0) - m2.get(k, 0.0) for k in keys]
    if norm == 'l1':
        return sum((abs(d) for d in diffs))
    return sqrt(sum((d * d for d in diffs)))

def quartet_distance(t1: Tree, t2: Tree) -> int:
    if t1.taxa != t2.taxa:
        raise ValueError('quartet_distance requires identical taxon sets (in order)')
    n = len(t1.taxa)
    ids = list(range(n))
    mismatches = 0
    for quad in combinations(ids, 4):
        if _quartet_topology(t1, quad) != _quartet_topology(t2, quad):
            mismatches += 1
    return mismatches

def _key(s: Split):
    return frozenset({s.A, s.B})

def _quartet_topology(t: Tree, quad):
    quad_set = frozenset(quad)
    for s in t.splits():
        a = s.A & quad_set
        b = s.B & quad_set
        if len(a) == 2 and len(b) == 2:
            return frozenset({frozenset(a), frozenset(b)})
    return 'star'