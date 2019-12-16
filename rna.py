from bisect import bisect_left
from graph import Node, build_components
from datagen import VCF
# from gene import make_genomic_graph, assign_reads_to_genomic_regions
from read import SNP, Read
from common import score
from rates import find_rates
from typing import Tuple, Dict, List, Set
from dataclasses import dataclass


@dataclass(init=False)
class Gene:
    id: int  # index
    name: str  # a.k.a. transcript
    chr: str
    interval: Tuple[int, int]
    sign: str

    exons: List[Tuple[int, int]]  # List of (start, len)
    neighbors: Set[int]
    snps: Set[int]

    def __init__(
        self,
        id: int,
        name: str,
        chr: str,
        interval: Tuple[int, int],
        sign: str,
        exons: List[Tuple[int, int]],
    ):
        self.id = id
        self.name = name
        self.chr = chr
        self.interval = interval
        self.sign = sign

        # Assumes that constructor calls sends in the start positions 1-based
        # (both for BED and GTF)
        self.exons = exons

        self.neighbors: Set[int] = set()
        self.snps: Set[int] = set()

    def __eq__(self, other):
        return self.name == other.name

    def set_snps(self, positions: List[Tuple[int, int]]):
        h = bisect_left(positions, self.interval[0])
        if h == len(positions):
            return
        exon = 0
        self.snps = set()
        for p in range(h, len(positions)):
            pos, snp = positions[p]
            if pos > self.interval[1]:
                break
            e_start, e_len = self.exons[exon]
            if pos >= e_start and pos < e_start + e_len - 1:
                self.snps.add(snp)
            elif pos >= start:
                exon += 1


@dataclass
class RNAGraph:
    genes: List[Gene]
    error: float
    ploidy: int
    single_reads: List[Read]
    multi_reads: List[Read]

    nodes: Dict[int, Node]


    comps: Dict[int, Set[int]]
    comps_dict: Dict[int, int]
    comps_mins_dict: Dict[int, Set[int]]
    SNP_to_genomic_region: Dict[int, int]
    genomic_region_to_SNPs: Dict[int, List[int]]
    reads_by_GR: Dict[int, List[Read]]
    # reads_by_indiv_gene: Dict[int, Dict[int, List[Read]]]
    # final: Dict[int, List[int]]
    snps_to_use: Dict[int, List[int]]

    reads_by_comsnps: Dict[int, List[Read]]
    read_dict: Dict[int, List[Read]]
    counts: Dict[int, List[int]]
    LLsnp: Dict[int, float]
    rates: Dict[int, Tuple[float, float]]

    components: Dict[int, List[int]]
    comp_mins: List[int]


    def __init__(
        self,
        vcf: VCF,
        genes: List[Gene],
        reads: List[Read],
        error: float
    ):
        self.genes = genes
        self.error = error
        self.ploidy = 2

        self.single_reads = [r for r in reads if len(r.snps) == 1]
        self.multi_reads = [r for r in reads if len(r.snps) > 1]

        chromosomes = {g.chr for g in genes}
        snps = sorted({k for r in self.single_reads for k in r.snps})
        self.nodes = {snp: Node(self.ploidy, snp, set()) for snp in snps}

        # Initialize gene SNPs and neighbours
        phasable_positions = {}  #S
        for s, node in self.nodes.items():
            snp = vcf.snps[s]
            phasable_positions.setdefault(snp.chr, []).append((snp.pos, snp))
        for p in phasable_positions.values():
            p.sort()
        for g in self.genes:
            if g.chr in phasable_positions:
                g.set_snps(phasable_positions[g.chr])
        s_dict: Dict[int, Set[int]] = {}  # SNP to gene
        for gi, g in enumerate(self.genes):
            for s in g.snps:
                s_dict[s].add(gi)
        for g in genes:
            g.neighbors = {n for s in g.snps for n in s_dict[s]}

        # Build components
        comps, comps_dict = build_components(genes)

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


        (
            self.comps,
            self.comps_dict,
            self.comps_mins_dict,
            self.SNP_to_genomic_region,
            self.genomic_region_to_SNPs,
            self.reads_by_GR
        ) = assign_reads_to_genomic_regions(self.genes, reads)
        # Here is where we decide which types of genes to use, those with no isoforms ("no_splicing")
        # all SNPs such that if there are multiple isoforms, those SNPs fall into all of them ("final")
        # self.final = self.dual_gene()
        self.snps_to_use = self.find_snps_to_use_all()

        self.reads_by_comsnps = self.find_reads_by_comsnps()
        self.read_dict = self.make_read_dict()  ##1-reads only
        self.counts = self.assign_counts_to_indiv_snps()
        self.LLsnp = self.assign_LL_dif_to_snps()
        print("assigning rates")
        self.rates = self.assign_rates2(self.snps_to_use)

    ######################################################################

    def find_snps_to_use_all(self) -> Dict[int, List[int]]:
        snps_to_use: Dict[int, List[int]] = {}
        for start in self.genomic_region_to_SNPs:
            snps = self.genomic_region_to_SNPs[start]
            if len(snps) > 1:
                snps_to_use[min(snps)] = sorted(snps)
        return snps_to_use

    def find_reads_by_comsnps(self) -> Dict[int, List[Read]]:
        # at this point we will only use reads that fall strictly within common snps
        # we should see how many long reads we arent using.
        # we should consider adding those back in to make components (earlier)
        reads_by_comsnps: Dict[int, List[Read]] = {}
        for start in sorted(self.snps_to_use.keys()):
            comp = self.snps_to_use[start]
            m = start
            reads_by_comsnps[start]: List[Read] = []
            gr = self.SNP_to_genomic_region[start]
            for r in self.reads_by_GR[gr]:
                tf = True
                for k in r.snps:
                    tf = tf and (k in comp)
                if tf:
                    reads_by_comsnps[start].append(r)
        return reads_by_comsnps

    def make_read_dict(self) -> Dict[int, List[Read]]:
        read_dicts: Dict[int, List[Read]] = {}
        for m in self.reads_by_comsnps:
            for r in self.reads_by_comsnps[m]:
                if len(r.snps) == 1:
                    k = r.snps[0]
                    if k in read_dicts:
                        read_dicts[k].append(r)
                    else:
                        read_dicts[k] = [r]
                else:
                    for k in r.snps:
                        new_read = Read({k: r.snps[k]}, r.count, -1)
                        if k in read_dicts:
                            read_dicts[k].append(new_read)
                        else:
                            read_dicts[k] = [new_read]
        return read_dicts

    def assign_counts_to_indiv_snps(self) -> Dict[int, List[int]]:
        ### ASSUMES 1READS ONLY #1reads
        all_counts: Dict[int, List[int]] = {}
        for s in self.read_dict:
            counts = [0, 0]
            for r in self.read_dict[s]:
                counts[r.snps[s] % 2] += r.count # originally r.read[s]
            all_counts[s] = counts
        return all_counts

    def assign_LL_dif_to_snps(self) -> Dict[int, float]:
        LL: Dict[int, float] = {}
        for start in self.snps_to_use:
            for snp in self.snps_to_use[start]:
                LL[snp] = score(self.counts[snp])
        return LL  # LL_gr,LL

    def assign_rates2(self, snps_to_use: Dict[int, List[int]]):
        rates: Dict[int, Tuple[float, float]] = {}
        for start in snps_to_use:
            snps = sorted(snps_to_use[start])
            reads = [r for s in snps for r in self.read_dict[s]]
            r = find_rates(snps, reads, 0.6)  ## Choose adjacent snp rate
            for s in snps:
                rates[s] = r
            for read in reads:
                read.rates = r
        return rates

    ############################################################

def make_data_from_RNA_data(RD: RNAData) -> Data:
    read_list = RD.multi_reads
    data = edges_from_readlist(read_list)
    return Data(
        data, RD.states, RD.k, RD.error, read_list, RD.positions, RD.names, RD.chroms
    )
