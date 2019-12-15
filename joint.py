from graph import Graph, mins_of_comps, build_all_comps
from read import Read
from rna import RNAData
from typing import Tuple, Dict, List, Set
from dataclasses import dataclass


@dataclass
class MiniNode:
    index: int
    neighbors: List[int]


@dataclass
class JointGraph:
    RD: RNAData
    G: Graph
    k: int
    short_comp_mins: List[int]
    short_comps_dict: Dict[int, int]
    long_comp_mins: List[int]
    nodes: Dict[int, MiniNode]
    nodekeys: List[int]

    components: Dict[int, List[int]]
    comp_mins: List[int]
    comps_dict: Dict[int, int]
    read_dict: Dict[int, List[Read]]
    comp_reads: Dict[int, List[Read]]

    forward_short: Dict[int, List[int]]
    forward_long: Dict[int, List[int]]
    back_short: Dict[int, int]
    back_long: Dict[int, int]


    def __init__(self, RD: RNAData, G: Graph):
        self.RD = RD
        self.G = G
        # self.read_list = RD.all_reads
        self.k = 2
        # self.error = G.error
        # self.short_components = RD.components  # components induced by genes  / DASE
        self.short_comp_mins = sorted(RD.comp_mins)
        self.short_comps_dict = {
            y: x for x in self.short_comp_mins for y in self.RD.components[x]
        }

        # components induced by 2+-reads (HapTree)
        self.long_comp_mins = sorted(G.comp_mins)
        # self.long_components = G.components
        # self.long_comps_dict = G.comps_dict

        # Organize nodes
        self.nodes: Dict[int, MiniNode] = {}
        Ls = len(self.short_comp_mins)
        Ll = len(self.long_comp_mins)

        self.forward_short: Dict[int, List[int]] = {}
        self.forward_long: Dict[int, List[int]] = {}
        self.back_short: Dict[int, int] = {}
        self.back_long: Dict[int, int] = {}
        for i in range(Ls):
            Ms = self.short_comp_mins[i]
            self.forward_short[i] = self.RD.components[Ms]
            self.back_short[Ms] = i
        for i in range(Ls, Ls + Ll):
            Ml = self.long_comp_mins[i - Ls]
            self.forward_long[i] = self.G.components[Ml]
            self.back_long[Ml] = i

        for i in range(Ls):
            self.nodes[i] = MiniNode(i, sorted(list({
                self.back_long[self.G.comps_dict[lil_node]]
                for lil_node in self.RD.components[self.short_comp_mins[i]]
                if lil_node in self.G.components
            })))

        for i in range(Ls, Ll + Ls):
            self.nodes[i] = MiniNode(i, sorted(list({
                self.back_short[self.short_comps_dict[lil_node]]
                for lil_node in self.G.components[self.long_comp_mins[i - Ls]]
                if lil_node in self.RD.components
            })))
        self.nodekeys = sorted(self.nodes.keys())

        # translate components
        weird_components, _ = build_all_comps(self)
        self.components: Dict[int, List[int]] = {}
        self.comps_dict: Dict[int, int] = {}
        for c in weird_components:
            temp_comp: Set[int] = {}
            for C in weird_components[c]:
                if C < Ls:
                    for x in self.forward_short[C]:
                        temp_comp.add(x)
                else:
                    for x in self.forward_long[C]:
                        temp_comp.add(x)
            comp = sorted(list(temp_comp))
            m = min(comp)
            for s in comp:
                self.components[s] = comp
                self.comps_dict[s] = m
        self.comp_mins = mins_of_comps(self.components)
        self.read_dict, self.comp_reads = self.make_read_dict()

    def make_read_dict(self) -> Tuple[Dict[int, List[Read]], Dict[int, List[Read]]]:
        # This does not work correctly if G is not the 2-reads from RNA-seq data.
        read_dict: Dict[int, List[Read]] = {}
        for s in self.RD.components:
            read_dict[s]: List[Read] = []
            for r in self.RD.read_dict[s]:  # 1-reads
                read_dict[s].append(r)
        comp_reads = {m: [] for m in self.comp_mins} #S
        for r in self.G.data.read_list.values():
            m = self.comps_dict[r.snps[0]]
            comp_reads[m].append(r)
            for k in r.snps:
                if k in read_dict:
                    read_dict[k].append(r)
                else:
                    read_dict[k] = [r]
        return read_dict, comp_reads

