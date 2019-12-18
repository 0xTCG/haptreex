from graph import Graph, Node, build_components, Component
from read import Read
from rna import RNAGraph
from typing import Tuple, Dict, List, Set
from dataclasses import dataclass
from pprint import pprint
import sys

@dataclass(init=False)
class JointGraph:
    ploidy: int
    components: Dict[int, Component]
    component_index: Dict[int, int]
    snp_reads: Dict[int, List[Read]]

    def __init__(self, rna: RNAGraph, dna: Graph):
        self.ploidy = rna.ploidy

        # Organize nodes
        forward_short, back_short = {}, {}
        for i, root in enumerate(rna.components):
            forward_short[i] = root
            back_short[root] = i
        forward_long, back_long = {}, {}
        for i, root in enumerate(dna.components):
            forward_long[i + len(rna.components)] = root
            back_long[root] = i + len(rna.components)
        nodes: Dict[int, Node] = {}
        for i, root in enumerate(rna.components):
            temp = {  # All DNA components that have some RNA nodes in this component
                back_long[dna.component_index[node]]
                for node in rna.components[root].nodes
                if node in dna.component_index
            }
            nodes[i] = Node(i, temp)
        for i, root in enumerate(dna.components):
            temp = {  # All DNA components that have some RNA nodes in this component
                back_short[rna.component_index[node]]
                for node in dna.components[root].nodes
                if node in rna.component_index
            }
            nodes[i + len(rna.components)] = Node(i + len(rna.components), temp)

        # Connect all DNA and RNA components...
        weird_components, _ = build_components(nodes)

        # translate components
        self.components: Dict[int, Component] = {}
        self.component_index: Dict[int, int] = {}
        for weird_root in weird_components:
            temp_comp: Set[int] = set()
            for weird_comp in weird_components[weird_root].nodes:
                if weird_comp < len(rna.components):
                    for x in rna.components[forward_short[weird_comp]].nodes:
                        temp_comp.add(x)
                else:
                    for x in dna.components[forward_long[weird_comp]].nodes:
                        temp_comp.add(x)
            n = sorted(temp_comp)
            root = n[0]
            self.components[root] = Component(root, n, [])
            for i in n:  # Update the dictionaries
                self.component_index[i] = root

        # This does not work correctly if G is not the 2-reads from RNA-seq data.
        self.snp_reads: Dict[int, List[Read]] = {}
        for snp in rna.component_index:  # 1-reads
            self.snp_reads[snp] = rna.snp_reads[snp][:]
        for read in dna.reads:
            for snp in read.snps:
                self.snp_reads.setdefault(snp, []).append(read)
