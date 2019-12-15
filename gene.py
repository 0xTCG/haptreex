from graph import Node
from read import Read
from typing import Tuple, Dict, List, Set
from dataclasses import dataclass


NOT_SEEN = -1
OVER_SEEN = -2


@dataclass
class Gene:
    name: str
    transcript_id: str
    chrom: str
    bp_start: int
    bp_end: int
    sign: str
    exon_count: int
    exon_lengths: List[int]
    exon_starts: List[int]
    gene_id: str
    gene_type: str

    ranges: List[Tuple[int, int]]
    neighbors: Set[int]
    snps: Set[int]
    index: int

    def __init__(
        self,
        name: str,
        chrom: str,
        bp_start: int,
        bp_end: int,
        sign: str,
        exon_count: int,
        exon_lengths: List[int],
        exon_starts: List[int],
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

        self.snps: Set[int] = {}
        self.neighbors: Set[int] = {}

        self.index = -1 # Dummy

    def __eq__(self, other):
        return self.name == other.name

    def __hash__(self):
        return self.name.__hash__()


def make_genomic_graph(GG: Dict[int, Gene]):  # making neighbours
    s_dict: Dict[int, Set[int]] = {} # snp to gene
    for g in GG.values():
        for s in g.snps:
            if s not in s_dict:
                s_dict[s] = {g.index}
            else:
                s_dict[s].add(g.index)
    for x in GG:
        GG[x].neighbors = {y for s in GG[x].snps for y in s_dict[s]}


def build_genomic_comp(start: Gene, GG: Dict[int, Gene]) -> Set[int]:
    # builds connected component containing the node start
    # returns ordered List of nodes in component
    s = {start.index}
    old_batch = {start.index}
    while len(old_batch) > 0:
        new_batch: Set[int] = {}
        for g in old_batch:
            new_batch = new_batch.union(GG[g].neighbors)
        new_batch = new_batch.difference(s)
        s = s.union(new_batch)
        old_batch = new_batch
    return s


def build_all_genomic_comps(
    GG: Dict[int, Gene]
) -> Tuple[Dict[int, Set[int]], Dict[int, int]]:
    ##For regular comps not genomic:
    # builds connected component dictionaries
    # comps maps each node to a List of nodes it is connected to
    # comps_dict maps each node to the smallest node it is connected to
    comps: Dict[int, Set[int]] = {}
    comps_dict: Dict[int, int] = {}
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
    comps: Dict[int, Set[int]], comps_dict: Dict[int, int]
) -> Dict[int, Set[int]]:
    mins = list(set(comps_dict.values()))
    gcbm: Dict[int, Set[int]] = {}
    for m in mins:
        gcbm[m] = comps[m]
    return gcbm


def snp_to_genomic_region(
    gcmb: Dict[int, Set[int]], GG: Dict[int, Gene]
) -> Tuple[Dict[int, int], Dict[int, List[int]]]:
    SNP_GR: Dict[int, int] = {}
    for m in gcmb:
        for g in gcmb[m]:
            for s in GG[g].snps:
                SNP_GR[s] = m
    GR_SNP = {m: [] for m in gcmb} #S
    for s in SNP_GR:
        GR_SNP[SNP_GR[s]].append(s)
    for m in GR_SNP:
        GR_SNP[m] = sorted(GR_SNP[m])
    return SNP_GR, GR_SNP


def assign_reads_to_genomic_regions(
    genes: Dict[int, Gene], reads: Dict[int, Read]
) -> Tuple[
    Dict[int, Set[int]],
    Dict[int, int],
    Dict[int, Set[int]],
    Dict[int, int],
    Dict[int, List[int]],
    Dict[int, List[Read]]
]:
    comps, comps_dict = build_all_genomic_comps(genes)
    comps_mins_dict = gene_comps_by_mins(comps, comps_dict)
    S_G, G_S = snp_to_genomic_region(comps_mins_dict, genes)
    reads_by_gene = {m: [] for m in comps_mins_dict} #S
    reads_by_gene[OVER_SEEN]: List[Read] = []
    reads_by_gene[NOT_SEEN]: List[Read] = []

    for R in reads.values():
        regions: Set[int] = {}
        for s in R.snps:
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


def determine_genes_gtf(gene_data: str, chroms: Set[str]) -> List[Gene]:
    print(f"Loading and formatting genes from {gene_data}...")

    b = [a[:-2].split("\t") for i, a in enumerate(open(gene_data, "r")) if a[0] != '#']
    trans_list: List[int] = []
    trans_dict: Dict[int, Dict[str, List[int]]] = {}

    for j, y in enumerate(b):
        if y[2] == "transcript":
            trans_list.append(j)
            trans_dict[j] = {"exon_starts": [], "exon_lengths": []} #S
    for k in range(len(trans_list) - 1):
        start = trans_list[k]
        end = trans_list[k + 1]
        for j in range(start, end):
            y = b[j]
            if y[2] == "exon":
                exon_start = int(y[3])
                exon_end = int(y[4])
                exon_length = exon_end - exon_start + 1
                trans_dict[start]["exon_starts"].append(exon_start)
                trans_dict[start]["exon_lengths"].append(exon_length)

    start = trans_list[-1]
    end = len(b)
    for j in range(start, end):
        y = b[j]
        if y[2] == "exon":
            exon_start = int(y[3])
            exon_end = int(y[4])
            exon_length = exon_end - exon_start + 1
            trans_dict[start]["exon_starts"].append(exon_start)
            trans_dict[start]["exon_lengths"].append(exon_length)

    transcripts: list[Gene] = {}
    for j in trans_dict:
        y = b[j]
        chrom = y[0]
        if chrom in chroms:
            bp_start = int(y[3])
            bp_end = int(y[4])
            sign = y[6]
            z = y[8]
            w = z.split("; ")
            v = [xx.split(" ") for xx in w]
            gene_id = v[0][1][1:-1]
            transcript_id = v[1][1][1:-1]
            gene_type = v[2][1][1:-1]
            exon_starts = trans_dict[j]["exon_starts"]
            exon_count = len(trans_dict[j]["exon_lengths"])
            exon_lengths = trans_dict[j]["exon_lengths"]
            transcripts[j].append(Gene(
                transcript_id,
                chrom,
                bp_start,
                bp_end,
                sign,
                exon_count,
                exon_lengths,
                exon_starts,
                gene_id,
                gene_type
            ))
            transcripts[-1].index = j
    return transcripts


def build_isodict(isoforms: str) -> Dict[str, Tuple[float, str]]:
    isodict: Dict[str, Tuple[float, str]] = {}
    with open(isoforms) as f:
        temp = [x.split("\t") for x in f.readlines()]
        for i in range(1, len(temp)):
            isodict[temp[i][0]] = float(temp[i][9]), temp[i][3]
    return isodict


def filter_transcripts(genes, isodict: Dict[str, Tuple[float, str]]) -> Dict[int, Gene]:
    gene_cov = {isodict[i][1]: 0.0 for i in isodict}
    for i in isodict:
        gene_cov[isodict[i][1]] += isodict[i][0]

    new_genes: Dict[int, Gene] = {}
    types: Set[str] = {}
    for j in genes:
        g = genes[j]
        tid = g.transcript_id
        types.add(genes[j].gene_type)
        if tid in isodict:
            val = isodict[tid][0]
            if val > 0:
                if val / gene_cov[isodict[tid][1]] > 0.01:
                    new_genes[j] = g
    return new_genes

