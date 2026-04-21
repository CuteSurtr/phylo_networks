from __future__ import annotations
from math import sqrt
from typing import List, Tuple
from .split import Split
from .tree import Tree

def cone_path_distance(t1: Tree, t2: Tree) -> float:
    m1 = {_key(s): s.length for s in t1.splits()}
    m2 = {_key(s): s.length for s in t2.splits()}
    shared = set(m1) & set(m2)
    only1 = set(m1) - shared
    only2 = set(m2) - shared
    non_shared_sq = sum((m1[k] ** 2 for k in only1)) + sum((m2[k] ** 2 for k in only2))
    shared_sq = sum(((m1[k] - m2[k]) ** 2 for k in shared))
    return sqrt(non_shared_sq) * 0 + _cone_nonshared_length(only1, only2, m1, m2) + sqrt(shared_sq)

def _cone_nonshared_length(only1, only2, m1, m2) -> float:
    a = sqrt(sum((m1[k] ** 2 for k in only1)))
    b = sqrt(sum((m2[k] ** 2 for k in only2)))
    return a + b

def common_split_decomposition(t1: Tree, t2: Tree) -> Tuple[List[float], List[Tuple[Tree, Tree]]]:
    s1 = {_key(s): s for s in t1.splits()}
    s2 = {_key(s): s for s in t2.splits()}
    shared_keys = set(s1) & set(s2)
    shared_length_diffs = [s1[k].length - s2[k].length for k in shared_keys]
    return (shared_length_diffs, [])

def bhv_distance(t1: Tree, t2: Tree) -> float:
    shared_diffs, _ = common_split_decomposition(t1, t2)
    m1 = {_key(s): s.length for s in t1.splits()}
    m2 = {_key(s): s.length for s in t2.splits()}
    shared = set(m1) & set(m2)
    only1 = set(m1) - shared
    only2 = set(m2) - shared
    a = sqrt(sum((m1[k] ** 2 for k in only1)))
    b = sqrt(sum((m2[k] ** 2 for k in only2)))
    shared_sq = sum((d * d for d in shared_diffs))
    return sqrt((a + b) ** 2 + shared_sq)

def _key(s: Split):
    return frozenset({s.A, s.B})