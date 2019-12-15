from itertools import combinations
from read import Read
from typing import Tuple, Dict, List, Set
from dataclasses import dataclass


def build_components(snps: List[SNP], nodes: Dict[SNP, Node]) -> Dict[SNP, List[SNP]]:
    """
    Build connected component dictionaries:
        {SNP: [sorted list of reachable SNPs]}.
    The component root is the first SNP.
    """
    comps: Dict[SNP, List[SNP]] = {}
    unbuilt = set(snps)
    while len(unbuilt) > 0:
        start = unbuilt.pop() # Take random SNP
        # Find all SNPs that are reachable from start
        comp, old_batch = {start}, {start}
        while len(old_batch) > 0:
            new_batch: Set[int] = {}
            for x in old_batch:
                new_batch |= set(nodes[x].neighbors)
            new_batch -= comp
            comp |= new_batch
            old_batch = new_batch
        unbuilt -= comp
        comp = sorted(comp)
        for snp in comp: # Update the dictionaries
            comps[snp] = comp
    return comps


@dataclass
class Node:
    ploidy: int
    snp: SNP
    neighbors: List[SNP]

    def __eq__(self, other):
        return self.snp == other.snp

    def __lt__(self, other):
        return self.snp < other.snp


@dataclass
class Edge:
    nodes: Tuple[Node, Node]
    ploidy: int

    def __eq__(self, other):
        return self.nodes == other.nodes


@dataclass
class Graph:
    reads: list[Read] # List of read objects
    ploidy: int
    error: float

    snps: List[SNP] # Sorted list of SNPs
    nodes: Dict[SNP, Node] # Mapping from a SNP to a node
    edges: Dict[Tuple[SNP, SNP], Edge] # TODO: even needed?

    component_roots: List[SNP]
    components: Dict[SNP, List[SNP]] # {SNP: [sorted list of reachable SNPs]}
    component_reads: Dict[SNP, List[Read]]
    snp_reads: Dict[SNP, List[Read]]

    def __init__(
        self,
        reads: List[Read],
        ploidy: int,
        error: float
    ):
        self.reads = reads
        self.ploidy = ploidy
        self.error = error

        data = { # All connected SNPs (edges)
            (i, j) if i < j else (j, i)
            for r in reads
            for i, j in combinations(r.snps, 2)
        }
        self.snps = sorted(list({i for e in data for i in e}))
        print(f"{len(self.snps)} SNPs in non-trivial connected components")
        # returns dictionary mapping NODE to those snps adjacent to NODE
        adj = {i: [] for i in self.snps} #S SNP -> adj. SNPs
        for e in data:
            adj[e[0]].append(e[1])
            adj[e[1]].append(e[0])
        self.nodes = {snp: Node(ploidy, snp, adj[snp]) for snp in self.snps}
        self.edges = {(i, j): Edge(self.ploidy, (self.nodes[i], self.nodes[j])) for i, j in data}

        """ReadGraph object, takes in data object"""

        self.components = build_components(self.snps, self.nodes)
        self.component_roots = [c[0] for c in self.components.values()]

        self.component_reads = {m: [] for m in self.comp_mins} #S Component root -> List of reads
        self.snp_reads = {s: [] for s in self.data.nodekeys} #S SNP -> List of reads
        for r in self.reads:
            self.component_reads.setdefault(self.components[r.snps[0]][0], []).append(r)
            for snp in r.snps[1:]:
                self.snp_reads.setdefault(snp, []).append(r)

    def __iter__(self):
        return self.nodes.values()
