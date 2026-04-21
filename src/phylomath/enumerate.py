from __future__ import annotations
from typing import Iterator, List
from .tree import Node, Tree, parse_newick, to_newick

def binary_topologies(n: int) -> Iterator[str]:
    if n < 1:
        return
    if n == 1:
        yield '1;'
        return
    trees: List[str] = ['(1,2);']
    for k in range(3, n + 1):
        next_trees: List[str] = []
        for s in trees:
            next_trees.extend(_insert_taxon(s, k))
        trees = next_trees
    yield from trees

def _insert_taxon(newick: str, k: int) -> Iterator[str]:
    body = newick[:-1] if newick.endswith(';') else newick
    n = len(body)
    starts: List[int] = []
    spans = []
    stack: List[int] = []
    i = 0
    while i < n:
        ch = body[i]
        if ch == '(':
            stack.append(i)
            i += 1
            continue
        if ch == ')':
            start = stack.pop()
            spans.append((start, i + 1))
            i += 1
            continue
        if ch == ',':
            i += 1
            continue
        j = i
        while j < n and body[j] not in ',)':
            j += 1
        spans.append((i, j))
        i = j
    seen = set()
    for start, end in spans:
        if (start, end) in seen:
            continue
        seen.add((start, end))
        sub = body[start:end]
        wrapped = body[:start] + f'({sub},{k})' + body[end:]
        yield (wrapped + ';')

def nni_neighbors(t: Tree) -> List[Tree]:
    results: List[Tree] = []
    base_newick = to_newick(t)
    base = parse_newick(base_newick)
    internal_edges: List[Node] = []
    for node in base.postorder():
        if node.parent is None or node.parent.parent is None:
            continue
        if node.is_leaf:
            continue
        internal_edges.append(node)
    for node in internal_edges:
        for swap_idx in (0, 1):
            trial = parse_newick(base_newick)
            _apply_nni(trial, _path_from_root(node), swap_idx)
            results.append(trial)
    return results

def _path_from_root(node: Node) -> List[int]:
    path: List[int] = []
    cur = node
    while cur.parent is not None:
        path.append(cur.parent.children.index(cur))
        cur = cur.parent
    return list(reversed(path))

def _follow(root: Node, path: List[int]) -> Node:
    n = root
    for i in path:
        n = n.children[i]
    return n

def _apply_nni(t: Tree, path: List[int], swap_idx: int) -> None:
    node = _follow(t.root, path)
    parent = node.parent
    if parent is None:
        return
    sibling = next((c for c in parent.children if c is not node))
    if len(node.children) < 2:
        return
    kid = node.children[swap_idx]
    p_idx = parent.children.index(sibling)
    n_idx = node.children.index(kid)
    parent.children[p_idx] = kid
    node.children[n_idx] = sibling
    kid.parent = parent
    sibling.parent = node