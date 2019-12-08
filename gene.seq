from graph import Node
from read import Read


NOT_SEEN = -1
OVER_SEEN = -2


class Gene:
    name: str
    transcript_id: str
    chrom: str
    bp_start: int
    bp_end: int
    sign: str
    exon_count: int
    exon_lengths: list[int]
    exon_starts: list[int]
    gene_id: str
    gene_type: str
    
    ranges: list[tuple[int, int]]
    neighbors: set[int]
    snps: set[int]
    index: int

    def __init__(
        self: Gene,
        name: str,
        chrom: str,
        bp_start: int,
        bp_end: int,
        sign: str,
        exon_count: int,
        exon_lengths: list[int],
        exon_starts: list[int],
        gene_id: str,
        gene_type: str
    ):
        self.name = name
        self.transcript_id = name
        self.gene_id = gene_id
        self.gene_type = gene_type
        self.chrom = chrom
        # Assumes that constructor calls sends in the start positions 1-based (both for BED and GTF)
        self.bp_start = bp_start
        self.bp_end = bp_end
        self.exon_lengths = exon_lengths
        self.exon_starts = exon_starts
        self.sign = sign
        self.exon_count = exon_count
        self.ranges = [
            (self.exon_starts[i], self.exon_starts[i] + self.exon_lengths[i] - 1)
            for i in range(exon_count)
        ]

        self.snps = set[int]()
        self.neighbors = set[int]()

        self.index = -1 # Dummy

    def __eq__(self: Gene, other: Gene):
        return self.name == other.name

    def __hash__(self: Gene):
        return self.name.__hash__()


def make_genomic_graph(GG: dict[int, Gene]):  # making neighbours
    s_dict = dict[int, set[int]]() # snp to gene
    for g in GG.values():
        for s in g.snps:
            if s not in s_dict:
                s_dict[s] = {g.index}
            else:
                s_dict[s].add(g.index)
    for x in GG:
        GG[x].neighbors = {y for s in GG[x].snps for y in s_dict[s]}


def build_genomic_comp(start: Gene, GG: dict[int, Gene]) -> set[int]:
    # builds connected component containing the node start
    # returns ordered list of nodes in component
    s = {start.index}
    old_batch = {start.index}
    while len(old_batch) > 0:
        new_batch = set[int]()
        for g in old_batch:
            new_batch = new_batch.union(GG[g].neighbors)
        new_batch = new_batch.difference(s)
        s = s.union(new_batch)
        old_batch = new_batch
    return s


def build_all_genomic_comps(
    GG: dict[int, Gene]
) -> tuple[dict[int, set[int]], dict[int, int]]:
    ##For regular comps not genomic:
    # builds connected component dictionaries
    # comps maps each node to a list of nodes it is connected to
    # comps_dict maps each node to the smallest node it is connected to
    comps = dict[int, set[int]]()
    comps_dict = dict[int, int]()
    for start in GG:
        if len(GG[start].snps) > 0:
            if start not in comps_dict:
                s = build_genomic_comp(GG[start], GG)
                comps[start] = s
                m = min(s)
                for g in s:
                    comps_dict[g] = m
                    comps[g] = s
    return comps, comps_dict


def gene_comps_by_mins(
    comps: dict[int, set[int]], comps_dict: dict[int, int]
) -> dict[int, set[int]]:
    mins = list(set(comps_dict.values()))
    gcbm = dict[int, set[int]]()
    for m in mins:
        gcbm[m] = comps[m]
    return gcbm


def snp_to_genomic_region(
    gcmb: dict[int, set[int]], GG: dict[int, Gene]
) -> tuple[dict[int, int], dict[int, list[int]]]:
    SNP_GR = dict[int, int]()
    for m in gcmb:
        for g in gcmb[m]:
            for s in GG[g].snps:
                SNP_GR[s] = m
    GR_SNP = {m: list[int]() for m in gcmb}
    for s in SNP_GR:
        GR_SNP[SNP_GR[s]].append(s)
    for m in GR_SNP:
        GR_SNP[m] = sorted(GR_SNP[m])
    return SNP_GR, GR_SNP


def assign_reads_to_genomic_regions(
    genes: dict[int, Gene], reads: dict[int, Read]
) -> tuple[
    dict[int, set[int]],
    dict[int, int],
    dict[int, set[int]],
    dict[int, int],
    dict[int, list[int]],
    dict[int, list[Read]]
]:
    comps, comps_dict = build_all_genomic_comps(genes)
    comps_mins_dict = gene_comps_by_mins(comps, comps_dict)
    S_G, G_S = snp_to_genomic_region(comps_mins_dict, genes)
    reads_by_gene = {m: list[Read]() for m in comps_mins_dict}
    reads_by_gene[OVER_SEEN] = list[Read]()
    reads_by_gene[NOT_SEEN] = list[Read]()

    for R in reads.values():
        regions = set[int]()
        for s in R.keys:
            if s in S_G:
                regions.add(S_G[s])
            else:
                regions.add(NOT_SEEN)
        if len(regions) == 1:
            m = list(regions)[0]
            # R.region = m
            reads_by_gene[m].append(R)
        else:
            reads_by_gene[OVER_SEEN].append(R)

    return comps, comps_dict, comps_mins_dict, S_G, G_S, reads_by_gene
