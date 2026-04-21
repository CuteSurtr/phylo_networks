from __future__ import annotations
from collections import deque
from typing import Dict, List
import networkx as nx
from .enumerate import binary_topologies, nni_neighbors
from .rf import rf_distance
from .tree import parse_newick, to_newick

def _canonical(newick: str) -> str:
    t = parse_newick(newick)

    def render(n) -> str:
        if n.is_leaf:
            return n.name or ''
        kid_strs = sorted((render(c) for c in n.children))
        return '(' + ','.join(kid_strs) + ')'
    return render(t.root) + ';'

def nni_graph(n: int) -> nx.Graph:
    G = nx.Graph()
    canonical_to_newick: Dict[str, str] = {}
    for s in binary_topologies(n):
        c = _canonical(s)
        canonical_to_newick[c] = s
    for c, s in canonical_to_newick.items():
        G.add_node(c)
    for c, s in canonical_to_newick.items():
        t = parse_newick(s)
        for nbr in nni_neighbors(t):
            nbr_canon = _canonical(to_newick(nbr))
            if nbr_canon in canonical_to_newick:
                G.add_edge(c, nbr_canon)
    return G

def rf_adjacency_bound_graph(n: int, max_rf: int=2) -> nx.Graph:
    G = nx.Graph()
    items = []
    for s in binary_topologies(n):
        c = _canonical(s)
        if c not in G:
            G.add_node(c)
            items.append((c, s))
    for i in range(len(items)):
        c1, s1 = items[i]
        t1 = parse_newick(s1)
        for j in range(i + 1, len(items)):
            c2, s2 = items[j]
            t2 = parse_newick(s2)
            if rf_distance(t1, t2) <= max_rf:
                G.add_edge(c1, c2)
    return G

def shortest_nni_distance(t1_newick: str, t2_newick: str, G: nx.Graph) -> int:
    c1, c2 = (_canonical(t1_newick), _canonical(t2_newick))
    return nx.shortest_path_length(G, c1, c2)