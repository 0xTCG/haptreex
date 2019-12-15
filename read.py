from random import choice
from typing import Tuple, Dict, List, Set, NamedTuple
from dataclasses import dataclass


class SNP(NamedTuple):
    chr: str
    pos: int
    name: str
    alleles: List[str]

    def __lt__(self, other):
        return (self.chr, self.pos) < (other.chr, other.pos)


@dataclass
class Read:
    id: int # Unique read ID
    count: int # Support
    snps: Dict[SNP, int] # SNP -> Allele
    # The following are used in alg.py:read_val_tail
    special_snp: SNP
    rates: List[float]

    def __init__(self, snps: Dict[SNP, int], count: int, id: int):
        self.id = id
        self.count = count
        self.snps = snps
        self.special_snp = min(self.snps)
        self.rates = [0.5, 0.5]

    def __eq__(self, other):
        return self.snps == other.snps and self.id == other.id

    def __len__(self):
        return len(self.snps)

    def __str__(self):
        return f'Read.{self.snps}'


def sample_from_reads(reads: List[Read]) -> List[Read]:
    return [
        r if len(r) else Read(dict(choice(r.snps.items())), 1, r.id)
        for r in reads
    ]
