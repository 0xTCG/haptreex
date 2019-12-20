from random import choice
from typing import Dict, List, NamedTuple
from dataclasses import dataclass


class SNP(NamedTuple):
    id: int
    chr: str
    pos: int
    name: str
    alleles: List[str]

    def __lt__(self, other):
        return self.id < other.id

    def __repr__(self):
        return f"{self.chr}:{self.pos + 1}"

    def __hash__(self):
        return id.__hash__()


@dataclass(init=False)
class Read:
    id: int  # Unique read ID
    count: int  # Support
    snps: Dict[int, int]  # SNP ID -> Allele
    # The following are used in alg.py:read_val_tail
    special_snp: int  # SNP ID
    rates: List[float]

    def __init__(self, snps: Dict[int, int], count: int, id: int):
        self.id = id
        self.count = count
        self.snps = snps
        self.special_snp = min(self.snps)
        self.rates = [0.5, 0.5]

    def __eq__(self, other):
        return self.id == other.id and self.snps == other.snps

    def __len__(self):
        return len(self.snps)

    def __str__(self):
        return f'Read.{self.snps}'


def sample_from_reads(reads: List[Read]) -> List[Read]:
    return [
        r if len(r) else Read(dict(choice(r.snps.items())), 1, r.id)
        for r in reads
    ]
