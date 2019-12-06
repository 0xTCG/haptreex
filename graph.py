from itertools import combinations
from read import Read

class Node:
    k: int
    index: int
    shape: int
    neighbors: list[int]
    pos: int
    name: str
    chrom: str
    num: int

    def __eq__(self: Node, other: Node):
        return self.index == other.index

    def __lt__(self: Node, other: Node):
        return self.index < other.index


class Edge:
    nodes: tuple[Node, Node]
    k: int
    indices: tuple[int, int]
    shapes: tuple[int, int]

    def __init__(self: Edge, node_pair: tuple[Node], k: int):
        self.nodes = node_pair
        self.k = k
        self.indices = (node_pair[0].index, node_pair[1].index)
        self.shapes = (node_pair[0].shape, node_pair[1].shape)

    def __eq__(self: Edge, other: Edge):
        return self.indices == other.indices


class Data:
    """
    A 'Data' object stores all initial info needed to create ReadGraph G
     - data contains read info
     - states is the state of each Node, in diploid this is always 1, counts the number of mutant alleles
     - k = ploidy
     - error = error rate assumed
     - read_list is list of reads
     - positions is the actual position on the genome (1-3bb) in base pair count
    """
    data: dict[tuple[int, int], int]
    states: dict[int, int]
    k: int
    error: float
    read_list: dict[int, Read]
    positions: dict[int, int]
    names: dict[int, str]
    chroms: dict[int, str]

    E: list[tuple[int, int]]
    nodekeys: list[int]
    ADJ: dict[int, list[int]]
    nodes: dict[int, Node]
    edges: dict[tuple[int, int], Edge]

    def __init__(
        self: Data,
        data: dict[tuple[int, int], int],
        states: dict[int, int],
        k: int,
        error: float,
        read_list: dict[int, Read],
        positions: dict[int, int],
        names: dict[int, str],
        chroms: dict[int, str],
    ):
        self.data = data
        self.states = states
        self.k = k
        self.error = error
        self.read_list = read_list
        self.positions = positions
        self.names = names
        self.chroms = chroms

        self.E = list(data.keys())
        self.nodekeys = self.node_keys()
        self.ADJ = self.adj()

        d = dict[int, Node]()
        for num in range(len(self.nodekeys)):
            index = self.nodekeys[num]
            d[index] = Node(
                self.k,
                index,
                self.states[index],
                self.ADJ[index],
                self.positions[index],
                self.names[index],
                self.chroms[index],
                num,
            )
        self.nodes = d

        N = self.nodes
        d = dict[tuple[int, int], Edge]()
        for e in self.E:
            d[e] = Edge([N[e[0]], N[e[1]]], self.k)
        self.edges = d

    def node_keys(self: Data):
        # computes all nodes actually seen by reads
        s = set[int]()
        for e in self.E:
            s.add(e[0])
            s.add(e[1])
        return sorted(list(s))

    def adj(self: Data):
        # returns dictionary mapping NODE to those nodes adjacent to NODE
        d = dict[int, list[int]]()
        for i in self.nodekeys:
            d[i] = list[int]()
        for e in self.E:
            d[e[0]].append(e[1])
            d[e[1]].append(e[0])
        return d


class Graph:
    data: Data
    size: int

    components: dict[int, list[int]]
    comps_dict: dict[int, int]
    comp_mins: list[int]

    POs: list[list[tuple[int, int]]]
    OPs: dict[tuple[int, int], int]
    read_dict: dict[int, list[Read]]
    comp_reads: dict[int, list[Read]]

    def __init__(self: Graph, data: Data):
        # ReadGraph object, takes in data object

        self.data = data
        self.size = len(self.data.nodekeys)
        print(f"{self.size} SNPs in non-trivial connected components")

        self.components, self.comps_dict = Graph.build_all_comps(self)
        self.comp_mins = mins_of_comps(self.components)

        self.POs = [[(0, 0)], [(0, 1), (1, 0)], [(1, 1)]]
        self.OPs = {(0, 1): 0, (1, 0): 1}
        self.read_dict, self.comp_reads = self.make_read_dict()

    def __iter__(self: Graph):
        return self.data.nodes.values()

    def make_read_dict(self) -> tuple[dict[int, list[Read]], dict[int, list[Read]]]:
        comp_reads = {m: list[Read]() for m in self.comp_mins}
        read_dicts = {s: list[Read]() for s in self.nodekeys}
        for r in self.read_list.values():
            m = self.comps_dict[r.keys[0]]
            comp_reads[m].append(r)
            for key in r.keys[1:]:
                read_dicts[key].append(r)
        return read_dicts, comp_reads

    def build_all_comps[G](g: G) -> tuple[dict[int, list[int]], dict[int, int]]:
        # builds connected component dictionaries
        # comps maps each node to a list of nodes it is connected to
        # comps_dict maps each node to the smallest node it is connected to
        comps = dict[int, list[int]]()
        comps_dict = dict[int, int]()
        unbuilt = set(g.nodekeys)
        while len(unbuilt) > 0:
            start = unbuilt.pop()
            s = build_comp(g, start)
            comps[start] = s
            comps_dict[start] = min(s)
            unbuilt.difference_update(s)  # we build_comp once for each set of components
            for key in s:  # and apply it to each nodekey
                comps[key] = s
                comps_dict[key] = comps_dict[start]
        return comps, comps_dict

    def build_comp[G](g: G, start: int) -> list[int]:
        # builds connected component containing the node start
        # returns ordered list of nodes in component
        s = set([start])
        old_batch = set([start])
        while len(old_batch) > 0:
            new_batch = set(list[int]())
            for x in old_batch:
                new_batch.update(set(g.nodes[x].neighbors))
            new_batch = new_batch.difference(s)
            s.update(new_batch)
            old_batch = new_batch

        return sorted(list(s))


def edges_from_readlist(read_list: dict[int, Read]):
    # returns dictionary with keys those edges that appear in G and values all 1
    d = dict[int, Read]()
    for zz in read_list:
        r = read_list[zz]
        seen_nodes = r.keys
        pairs = list(combinations(seen_nodes, 2))
        for pair in pairs:
            i, j = pair
            pair = (i, j) if i < j else (j, i)
            d[pair] = 1
    return d


def mins_of_comps(comps: dict[int, list[int]]) -> list[int]:
    mins = dict[int, int]()
    for comp in comps.values():
        mins[min(comp)] = 1
    return sorted(mins.keys())

