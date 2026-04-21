from __future__ import annotations
from dataclasses import dataclass, field
from typing import FrozenSet, Iterable

@dataclass(frozen=True)
class Split:
    A: FrozenSet[int]
    B: FrozenSet[int]
    length: float = field(default=0.0, compare=False, hash=False)

    @staticmethod
    def of(block: Iterable[int], taxa: Iterable[int], length: float=0.0) -> 'Split':
        block_set = frozenset(block)
        taxa_set = frozenset(taxa)
        complement = taxa_set - block_set
        if not block_set or not complement:
            raise ValueError('split blocks must both be non-empty')
        ref = min(taxa_set)
        if ref in block_set:
            A, B = (complement, block_set)
        else:
            A, B = (block_set, complement)
        return Split(A=A, B=B, length=float(length))

    @property
    def is_trivial(self) -> bool:
        return len(self.A) == 1 or len(self.B) == 1

    def compatible_with(self, other: 'Split') -> bool:
        if not self.A & other.A:
            return True
        if not self.A & other.B:
            return True
        if not self.B & other.A:
            return True
        if not self.B & other.B:
            return True
        return False

    def __repr__(self) -> str:
        a = '{' + ','.join(map(str, sorted(self.A))) + '}'
        b = '{' + ','.join(map(str, sorted(self.B))) + '}'
        return f'Split({a}|{b}, ℓ={self.length:g})'

def pairwise_compatible(splits: Iterable[Split]) -> bool:
    sp = list(splits)
    for i in range(len(sp)):
        for j in range(i + 1, len(sp)):
            if not sp[i].compatible_with(sp[j]):
                return False
    return True