from graph import Graph, mins_of_comps
from read import Read
from rna import RNAData


class MiniNode:
    index: int
    neighbors: list[int]


class JointGraph:
    RD: RNAData
    G: Graph
    k: int
    short_comp_mins: list[int]
    short_comps_dict: dict[int, int]
    long_comp_mins: list[int]
    nodes: dict[int, MiniNode]
    nodekeys: list[int]

    weird_components: dict[int, list[int]]
    components: dict[int, list[int]]
    comp_mins: list[int]()
    comps_dict: dict[int, int]
    read_dict: dict[int, list[Read]]
    comp_reads: dict[int, list[Read]]

    def __init__(self: JointGraph, RD: RNAData, G: Graph):
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

        self.nodes = self.organize_nodes()
        self.nodekeys = sorted(self.nodes.keys())

        self.weird_components, _ = Graph.build_all_comps(self)

        self.components, self.comp_mins, self.comps_dict = self.translate_components()
        self.comp_mins = mins_of_comps(self.components)
        self.read_dict, self.comp_reads = self.make_read_dict()

    def organize_nodes(self: JointGraph) -> dict[int, MiniNode]:
        nodes = {}
        Ls = len(self.short_comp_mins)
        self.Ls = Ls
        Ll = len(self.long_comp_mins)
        self.Ll = Ll

        forward_short = {}
        forward_long = {}
        back_short = {}
        back_long = {}
        for i in range(Ls):
            Ms = self.short_comp_mins[i]
            forward_short[i] = self.RD.components[Ms]
            back_short[Ms] = i
        for i in range(Ls, Ls + Ll):
            Ml = self.long_comp_mins[i - Ls]
            forward_long[i] = self.long_components[Ml]
            back_long[Ml] = i

        for i in range(Ls):
            temp_neighbors = set()
            for lil_node in self.RD.components[self.short_comp_mins[i]]:
                if lil_node in self.long_components:
                    temp_neighbors.add(back_long[self.long_comps_dict[lil_node]])

            nodes[i] = mini_node(i, sorted(list(temp_neighbors)))

        for i in range(Ls, Ll + Ls):
            temp_neighbors = set()
            for lil_node in self.long_components[self.long_comp_mins[i - Ls]]:
                if lil_node in self.RD.components:
                    temp_neighbors.add(back_short[self.short_comps_dict[lil_node]])

            nodes[i] = mini_node(i, sorted(list(temp_neighbors)))

        self.forward_short = forward_short
        self.forward_long = forward_long
        self.back_short = back_short
        self.back_long = back_long

        return nodes

    def translate_components(
        self: JointGraph,
    ) -> tuple[dict[int, list[int]], list[int], dict[int, int]]:
        comp_mins = list[int]()
        components = dict[int, list[int]]()
        comps_dict = dict[int, int]()
        for c in self.weird_components:
            temp_comp = set()
            for C in self.weird_components[c]:
                if C < self.Ls:
                    for x in self.forward_short[C]:
                        temp_comp.add(x)
                else:
                    for x in self.forward_long[C]:
                        temp_comp.add(x)
            comp = sorted(list(temp_comp))
            m = min(comp)
            comp_mins.append(m)
            for s in comp:
                components[s] = comp
                comps_dict[s] = m
        return components, sorted(comp_mins), comps_dict

    def make_read_dict(
        self: JointGraph
    ) -> tuple[dict[int, list[Read]], dict[int, list[Read]]]:
        # This does not work correctly if G is not the 2-reads from RNA-seq data.
        read_dict = dict[int, list[Read]]()
        for s in self.RD.components:
            read_dict[s] = list[Read]()
            for r in self.RD.read_dict[s]:  # 1-reads
                read_dict[s].append(r)
        comp_reads = {m: list[Read]() for m in self.comp_mins}
        for r in self.G.read_list.values():
            m = self.comps_dict[r.keys[0]]
            comp_reads[m].append(r)
            for k in r.keys:
                if k in read_dict:
                    read_dict[k].append(r)
                else:
                    read_dict[k] = [r]
        return read_dict, comp_reads
