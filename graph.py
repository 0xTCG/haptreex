from itertools import combinations
from read import SNP, Read
from typing import Tuple, Dict, List, Set
from dataclasses import dataclass


@dataclass
class Node:
    ploidy: int
    snp: int  # Here, only SNP IDs are being used
    neighbors: Set[int]

    def __eq__(self, other):
        return self.snp == other.snp

    def __lt__(self, other):
        return self.snp < other.snp


@dataclass
class Component:
    root: int  # Root SNP
    snps: List[int]  # List of nodes
    reads: List[Read]  # Supporting reads

    # component_roots: List[int]  # List of root SNP IDs for each component
    # components: Dict[int, List[int]]  # {SNP ID: [sorted list of reachable SNP IDs]}
    # component_reads: Dict[int, List[Read]]  # List of reads that support component
    # snp_reads: Dict[int, List[Read]]  # Reads for each SNP



def build_components(
    snps: List[int],
    nodes: Dict[int, Node]
) -> Tuple[Dict[int, Component], Dict[int, int]]:
    """
    Build connected components:
        {SNP ID: Component}
    and the associated index:
        {SNP ID: component root ID}.
    The root is the smallest SNP in the component.
    """
    comps: Dict[int, Component] = {}
    index: Dict[int, int] = {}
    unbuilt = set(snps)
    while len(unbuilt) > 0:
        start = unbuilt.pop()  # Take a random SNP
        # Find all SNPs that are reachable from start
        comp, old_batch = {start}, {start}
        while len(old_batch) > 0:
            new_batch: Set[int] = set()
            for x in old_batch:
                new_batch |= nodes[x].neighbors
            new_batch -= comp
            comp |= new_batch
            old_batch = new_batch
        unbuilt -= comp

        snps = sorted(comp)
        root = snps[0]
        comps[root] = Component(root, snps, [])
        for snp in snps:  # Update the dictionaries
            index[snp] = root
    return comps, index


@dataclass(init=False)
class Graph:
    reads: List[Read]  # List of read objects
    ploidy: int
    error: float

    snps: List[int]  # Sorted list of SNP IDs
    nodes: Dict[int, Node]  # Mapping from a SNP ID to a node

    components: Dict[int, Component]  # {Component ID: Component}
    component_index: Dict[int, int]  # SNP ID: Component ID
    snp_reads: Dict[int, List[Read]]  # Reads for each SNP

    def __init__(
        self,
        reads: List[Read],
        ploidy: int,
        error: float
    ):
        self.reads = reads
        self.ploidy = ploidy
        self.error = error

        self.snps = sorted({s for r in reads for s in r.snps})
        print(f"{len(self.snps)} SNPs in non-trivial connected components")
        # returns dictionary mapping NODE to those snps adjacent to NODE
        adj: Dict[int, Set[int]] = {i: set() for i in self.snps}  #S SNP -> adj. SNPs
        for r in reads:
            for i in r.snps:
                for j in r.snps:
                    adj[i].add(j)
        self.nodes = {snp: Node(ploidy, snp, adj[snp]) for snp in self.snps}

        # Build connected components
        self.components, self.component_index = build_components(self.snps, self.nodes)

        self.snp_reads = {s: [] for s in self.snps}  #S SNP -> List of reads
        for r in self.reads:
            root = self.component_index[r.special_snp]
            self.components[root].reads.append(r)
            for snp in r.snps:
                if snp != root:
                    self.snp_reads.setdefault(snp, []).append(r)

    def __iter__(self):
        return self.nodes.values()
