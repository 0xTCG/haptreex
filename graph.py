from itertools import combinations
from read import SNP, Read
from typing import Tuple, Dict, List, Set, Any
from dataclasses import dataclass
# from rna import RNAGraph


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

    def integrate_rna(self, rna):
        assert self.ploidy == 2

        # Organize nodes
        forward_short, back_short = {}, {}
        for i, root in enumerate(rna.components):
            forward_short[i] = root
            back_short[root] = i
        forward_long, back_long = {}, {}
        for i, root in enumerate(self.components):
            forward_long[i + len(rna.components)] = root
            back_long[root] = i + len(rna.components)
        nodes: Dict[int, Node] = {}
        for i, root in enumerate(rna.components):
            temp = {  # All DNA components that have some RNA nodes in this component
                back_long[self.component_index[node]]
                for node in rna.components[root].nodes
                if node in self.component_index
            }
            nodes[i] = Node(i, temp)
        for i, root in enumerate(self.components):
            temp = {  # All DNA components that have some RNA nodes in this component
                back_short[rna.component_index[node]]
                for node in self.components[root].nodes
                if node in rna.component_index
            }
            nodes[i + len(rna.components)] = Node(i + len(rna.components), temp)

        # Connect all DNA and RNA components...
        weird_components, _ = build_components(nodes)

        # translate components
        components: Dict[int, Component] = {}
        component_index: Dict[int, int] = {}
        for weird_root in weird_components:
            temp_comp: Set[int] = set()
            for weird_comp in weird_components[weird_root].nodes:
                if weird_comp < len(rna.components):
                    for x in rna.components[forward_short[weird_comp]].nodes:
                        temp_comp.add(x)
                else:
                    for x in self.components[forward_long[weird_comp]].nodes:
                        temp_comp.add(x)
            n = sorted(temp_comp)
            root = n[0]
            components[root] = Component(root, n, [])
            for i in n:  # Update the dictionaries
                component_index[i] = root
        self.components, self.component_index = components, component_index

        # This does not work correctly if G is not the 2-reads from RNA-seq data.
        self.snp_reads: Dict[int, List[Read]] = {}
        for snp in rna.component_index:  # 1-reads
            self.snp_reads[snp] = rna.snp_reads[snp][:]
        for read in self.reads:
            for snp in read.snps:
                self.snp_reads.setdefault(snp, []).append(read)
