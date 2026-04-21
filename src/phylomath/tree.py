from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, FrozenSet, List, Optional, Set, Tuple
from .split import Split

@dataclass
class Node:
    name: Optional[str] = None
    taxon_id: Optional[int] = None
    length: float = 0.0
    parent: Optional['Node'] = None
    children: List['Node'] = field(default_factory=list)

    @property
    def is_leaf(self) -> bool:
        return not self.children

@dataclass
class Tree:
    root: Node
    taxa: List[str]
    _label_to_id: Dict[str, int] = field(default_factory=dict)

    @staticmethod
    def build(root: Node, taxa: List[str]) -> 'Tree':
        canonical = sorted(set(taxa))
        lookup = {t: i for i, t in enumerate(canonical)}
        t = Tree(root=root, taxa=canonical, _label_to_id=lookup)
        t._assign_taxon_ids()
        return t

    def _assign_taxon_ids(self) -> None:
        for n in self.postorder():
            if n.is_leaf and n.name is not None:
                if n.name not in self._label_to_id:
                    raise ValueError(f'unknown taxon label {n.name!r}')
                n.taxon_id = self._label_to_id[n.name]

    def postorder(self) -> List[Node]:
        out: List[Node] = []

        def visit(n: Node) -> None:
            for c in n.children:
                visit(c)
            out.append(n)
        visit(self.root)
        return out

    def leaves(self) -> List[Node]:
        return [n for n in self.postorder() if n.is_leaf]

    def splits(self, include_trivial: bool=False) -> List[Split]:
        taxa_ids = frozenset(range(len(self.taxa)))
        below: Dict[int, FrozenSet[int]] = {}
        from collections import OrderedDict
        length_by_key: 'OrderedDict[Tuple, Tuple[Split, float]]' = OrderedDict()
        for n in self.postorder():
            if n.is_leaf:
                below[id(n)] = frozenset({n.taxon_id})
                continue
            acc: Set[int] = set()
            for c in n.children:
                acc |= below[id(c)]
            below[id(n)] = frozenset(acc)
            if n.parent is None:
                continue
            s = Split.of(acc, taxa_ids, length=n.length)
            if not include_trivial and s.is_trivial:
                continue
            key = frozenset({s.A, s.B})
            if key in length_by_key:
                prev, acc_len = length_by_key[key]
                length_by_key[key] = (prev, acc_len + n.length)
            else:
                length_by_key[key] = (s, n.length)
        merged: List[Split] = []
        for _, (s, total_len) in length_by_key.items():
            merged.append(Split(A=s.A, B=s.B, length=total_len))
        if include_trivial:
            for lf in self.leaves():
                s = Split.of({lf.taxon_id}, taxa_ids, length=lf.length)
                merged.append(s)
        return merged

def parse_newick(s: str) -> Tree:
    if not isinstance(s, str):
        raise TypeError(f'parse_newick expected a str, got {type(s).__name__}')
    s_stripped = s.strip()
    if not s_stripped:
        raise ValueError('parse_newick: empty input')
    if '(' not in s_stripped and ';' not in s_stripped:
        raise ValueError(f'parse_newick: input does not look like Newick: {s_stripped!r}')
    tokens = _tokenize(s_stripped)
    pos = [0]

    def peek() -> str:
        return tokens[pos[0]]

    def eat(tok: str) -> None:
        if tokens[pos[0]] != tok:
            raise ValueError(f'expected {tok!r} at position {pos[0]}, got {tokens[pos[0]]!r}')
        pos[0] += 1

    def parse_subtree() -> Node:
        if peek() == '(':
            eat('(')
            children = [parse_subtree()]
            while peek() == ',':
                eat(',')
                children.append(parse_subtree())
            eat(')')
            name = None
            if peek() not in (',', ')', ':', ';'):
                name = tokens[pos[0]]
                pos[0] += 1
            length = 0.0
            if peek() == ':':
                eat(':')
                length = float(tokens[pos[0]])
                pos[0] += 1
            node = Node(name=name, length=length, children=children)
            for c in children:
                c.parent = node
            return node
        name = tokens[pos[0]]
        pos[0] += 1
        length = 0.0
        if peek() == ':':
            eat(':')
            length = float(tokens[pos[0]])
            pos[0] += 1
        return Node(name=name, length=length)
    root = parse_subtree()
    if peek() == ';':
        eat(';')
    taxa: List[str] = []
    seen: Set[str] = set()

    def collect(n: Node) -> None:
        if n.is_leaf and n.name is not None and (n.name not in seen):
            seen.add(n.name)
            taxa.append(n.name)
        for c in n.children:
            collect(c)
    collect(root)
    return Tree.build(root, taxa)

def _tokenize(s: str) -> List[str]:
    out: List[str] = []
    i, n = (0, len(s))
    while i < n:
        c = s[i]
        if c.isspace():
            i += 1
            continue
        if c in '(),:;':
            out.append(c)
            i += 1
            continue
        j = i
        while j < n and s[j] not in '(),:;' and (not s[j].isspace()):
            j += 1
        out.append(s[i:j])
        i = j
    out.append(';' if out and out[-1] != ';' else ';')
    return out

def to_newick(t: Tree) -> str:

    def render(n: Node) -> str:
        if n.is_leaf:
            body = n.name or ''
        else:
            body = '(' + ','.join((render(c) for c in n.children)) + ')'
            if n.name:
                body += n.name
        if n.parent is not None and n.length != 0.0:
            body += f':{n.length:g}'
        return body
    return render(t.root) + ';'