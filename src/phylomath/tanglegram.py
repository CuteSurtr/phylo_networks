from __future__ import annotations
from itertools import product
from typing import Dict, List, Optional, Tuple
from .tree import Node, Tree

def leaf_order(t: Tree, flip_mask: Optional[List[bool]]=None) -> List[int]:
    if flip_mask is None:
        return [lf.taxon_id for lf in _inorder_leaves(t.root)]
    internal: List[Node] = [n for n in t.postorder() if not n.is_leaf]
    flips = dict(zip((id(n) for n in internal), flip_mask))

    def walk(n: Node) -> List[int]:
        if n.is_leaf:
            return [n.taxon_id]
        kids = list(n.children)
        if flips.get(id(n), False):
            kids = list(reversed(kids))
        out: List[int] = []
        for c in kids:
            out.extend(walk(c))
        return out
    return walk(t.root)

def _inorder_leaves(n: Node) -> List[Node]:
    if n.is_leaf:
        return [n]
    out: List[Node] = []
    for c in n.children:
        out.extend(_inorder_leaves(c))
    return out

def count_crossings(order1: List[int], order2: List[int]) -> int:
    if sorted(order1) != sorted(order2):
        raise ValueError('leaf orderings must be permutations of the same taxon set')
    pos_in_1 = {tid: i for i, tid in enumerate(order1)}
    perm = [pos_in_1[tid] for tid in order2]
    return _count_inversions(perm)

def _count_inversions(a: List[int]) -> int:

    def merge_count(arr: List[int]) -> Tuple[List[int], int]:
        if len(arr) <= 1:
            return (arr, 0)
        mid = len(arr) // 2
        left, lc = merge_count(arr[:mid])
        right, rc = merge_count(arr[mid:])
        merged: List[int] = []
        i = j = inv = 0
        while i < len(left) and j < len(right):
            if left[i] <= right[j]:
                merged.append(left[i])
                i += 1
            else:
                merged.append(right[j])
                inv += len(left) - i
                j += 1
        merged.extend(left[i:])
        merged.extend(right[j:])
        return (merged, lc + rc + inv)
    _, inv = merge_count(list(a))
    return inv

def min_crossings_one_tree(t_fixed: Tree, t_free: Tree) -> int:
    if sorted(t_fixed.taxa) != sorted(t_free.taxa):
        raise ValueError('trees must have the same taxon set')
    pos_in_fixed = {tid: i for i, tid in enumerate(leaf_order(t_fixed))}

    def rec(n) -> tuple:
        if n.is_leaf:
            return (0, [n.taxon_id])
        child_results = [rec(c) for c in n.children]
        if len(child_results) == 2:
            cL, oL = child_results[0]
            cR, oR = child_results[1]
            cost_noflip = _cross_between(oL, oR, pos_in_fixed)
            cost_flip = _cross_between(oR, oL, pos_in_fixed)
            if cost_noflip <= cost_flip:
                return (cL + cR + cost_noflip, oL + oR)
            return (cL + cR + cost_flip, oR + oL)
        best = None
        best_order = None
        from itertools import permutations
        for perm in permutations(range(len(child_results))):
            order: list = []
            cost = sum((child_results[i][0] for i in perm))
            for i in range(len(perm)):
                for j in range(i + 1, len(perm)):
                    cost += _cross_between(child_results[perm[i]][1], child_results[perm[j]][1], pos_in_fixed)
                order += child_results[perm[i]][1]
            if best is None or cost < best:
                best = cost
                best_order = order + child_results[perm[-1]][1]
        return (best or 0, best_order or [])
    total, _ = rec(t_free.root)
    return total

def _cross_between(left_tids, right_tids, pos_in_fixed: dict) -> int:
    crossings = 0
    for a in left_tids:
        pa = pos_in_fixed[a]
        for b in right_tids:
            if pa > pos_in_fixed[b]:
                crossings += 1
    return crossings

def min_crossings_exhaustive(t1: Tree, t2: Tree) -> int:
    if sorted(t1.taxa) != sorted(t2.taxa):
        raise ValueError('trees must have the same taxon set')
    in1 = [n for n in t1.postorder() if not n.is_leaf]
    in2 = [n for n in t2.postorder() if not n.is_leaf]
    best = None
    for m1 in product((False, True), repeat=len(in1)):
        o1 = leaf_order(t1, list(m1))
        for m2 in product((False, True), repeat=len(in2)):
            o2 = leaf_order(t2, list(m2))
            c = count_crossings(o1, o2)
            if best is None or c < best:
                best = c
    return best if best is not None else 0