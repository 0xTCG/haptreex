from itertools import combinations
from read import SNP, Read
from typing import Tuple, Dict, List, Set, Any
from dataclasses import dataclass


@dataclass
class Node:
    snp: int  # Here, only SNP IDs are being used
    neighbors: Set[int]

    def __eq__(self, other):
        return self.snp == other.snp

    def __lt__(self, other):
        return self.snp < other.snp


@dataclass
class Component:
    root: int  # Root node ID
    nodes: List[int]  # List of nodes
    reads: List[Read]  # Supporting reads


def build_components(
    nodes: Dict[int, Any]
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
    unbuilt = set(nodes.keys())
    while len(unbuilt) > 0:
        start = unbuilt.pop()  # Take a random SNP
        # Find all nodes that are reachable from start
        comp, old_batch = {start}, {start}
        while len(old_batch) > 0:
            new_batch: Set[int] = set()
            for x in old_batch:
                new_batch |= nodes[x].neighbors
            new_batch -= comp
            comp |= new_batch
            old_batch = new_batch
        unbuilt -= comp

        n = sorted(comp)
        root = n[0]
        comps[root] = Component(root, n, [])
        for i in n:  # Update the dictionaries
            index[i] = root
    return comps, index


@dataclass(init=False)
class Graph:
    reads: List[Read]  # List of read objects
    ploidy: int

    snps: List[int]  # Sorted list of SNP IDs
    nodes: Dict[int, Node]  # Mapping from a SNP ID to a node

    components: Dict[int, Component]  # {Component ID: Component}
    component_index: Dict[int, int]
    snp_reads: Dict[int, List[Read]]  # Reads for each SNP

    def __init__(
        self,
        reads: List[Read],
        ploidy: int
    ):
        self.reads = reads
        self.ploidy = ploidy

        self.snps = sorted({s for r in reads for s in r.snps})
        print(f"{len(self.snps)} SNPs in non-trivial connected components")
        adj: Dict[int, Set[int]] = {i: set() for i in self.snps}  #S SNP -> adj. SNPs
        for r in reads:
            for i in r.snps:
                for j in r.snps:
                    adj[i].add(j)
        self.nodes = {snp: Node(snp, adj[snp]) for snp in self.snps}

        # Build connected components
        self.components, self.component_index = build_components(self.nodes)

        self.snp_reads = {s: [] for s in self.snps}  #S SNP -> List of reads
        for r in self.reads:
            root = self.component_index[r.special_snp]
            self.components[root].reads.append(r)
            for snp in r.snps:
                if snp != root:
                    self.snp_reads.setdefault(snp, []).append(r)

    def __iter__(self):
        return self.nodes.values()
