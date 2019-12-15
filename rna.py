from bisect import bisect_left
from graph import Node, Data, edges_from_readlist
from gene import Gene, make_genomic_graph, assign_reads_to_genomic_regions
from read import Read
from common import score
from rates import find_rates
from typing import Tuple, Dict, List, Set
from dataclasses import dataclass


@dataclass
class RNAData:
    states: Dict[int, int]
    # initial_genes: Dict[int, Gene]
    genes: Dict[int, Gene]
    error: float
    # all_reads: Dict[int, Read]
    positions: Dict[int, int]
    names: Dict[int, str]
    chroms: Dict[int, str]
    isodict: Dict[str, Tuple[float, str]]
    k: int

    read_list: Dict[int, Read]
    long_read_list: Dict[int, Read]

    nodekeys: List[int]
    nodes: Dict[int, Node]
    chrom_list: List[str]
    phasable_positions: Dict[str, List[int]]
    PSdict: Dict[str, Dict[int, int]]

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
        states: Dict[int, int],
        genes: Dict[int, Gene],
        filtered_genes: Dict[int, Gene],
        error: float,
        read_list: Dict[int, Read],
        positions: Dict[int, int],
        names: Dict[int, str],
        chroms: Dict[int, str],
        isodict: Dict[str, Tuple[float, str]]
    ):
        self.k = 2
        self.states = states
        # self.initial_genes = genes
        self.genes = filtered_genes
        self.error = error
        # self.all_reads = read_list
        self.positions = positions
        self.names = names
        self.chroms = chroms
        self.isodict = isodict

        short_read_list: Dict[int, Read] = {}
        long_read_list: Dict[int, Read] = {}
        for i in read_list:
            r = read_list[i]
            r.rates = [0.5, 0.5]
            if len(r.snps) == 1:
                short_read_list[i] = r
            else:
                long_read_list[i] = r

        self.read_list = short_read_list
        # self.short_read_list = short_read_list
        self.long_read_list = long_read_list

        self.nodekeys = sorted(list(self.node_keys()))

        d: Dict[int, Node] = {}
        chrom_list: Dict[str, int] = {}
        for num in range(len(self.nodekeys)):
            index = self.nodekeys[num]
            d[index] = Node(
                self.k,
                index,
                0,
                [], #S List[int](),
                self.positions[index],
                self.names[index],
                self.chroms[index],
                num
            )
            chrom_list[self.chroms[index]] = 0
        self.nodes = d
        self.chrom_list = sorted(chrom_list.keys())
        self.phasable_positions = {chrom: [] for chrom in self.chrom_list} #S
        for x in self.nodekeys:
            self.phasable_positions[self.chroms[x]].append(self.positions[x])
        self.PSdict = self.make_position_dict()

        ######################################################################

        for g in self.genes:
            self.find_SNPs(self.genes[g])
        make_genomic_graph(self.genes)

        ######################################################################

        reads = read_list
        (
            self.comps,
            self.comps_dict,
            self.comps_mins_dict,
            self.SNP_to_genomic_region,
            self.genomic_region_to_SNPs,
            self.reads_by_GR
        ) = assign_reads_to_genomic_regions(self.genes, reads)
        # self.reads_by_indiv_gene = self.assign_reads_to_indiv_genes()

        ######################################################################

        # Here is where we decide which types of genes to use, those with no isoforms ("no_splicing")
        # all SNPs such that if there are multiple isoforms, those SNPs fall into all of them ("final")

        # self.final = self.dual_gene()
        self.snps_to_use = self.find_snps_to_use_all()

        ######################################################################

        self.reads_by_comsnps = self.find_reads_by_comsnps()
        self.read_dict = self.make_read_dict()  ##1-reads only
        self.counts = self.assign_counts_to_indiv_snps()
        self.LLsnp = self.assign_LL_dif_to_snps()
        print("assigning rates")
        self.rates = self.assign_rates2(self.snps_to_use)

    ######################################################################

    def node_keys(self) -> List[int]:
        s = {k: 0 for r in self.read_list.values() for k in r.snps}
        return list(s.keys())

    def make_position_dict(self) -> Dict[str, Dict[int, int]]:
        position_to_SNP_index_dict = {chrom: Dict[int, int]() for chrom in self.chrom_list}
        for s in self.nodekeys:
            S = self.nodes[s]
            c = S.chrom
            p = S.pos
            position_to_SNP_index_dict[c][p] = S.index
        return position_to_SNP_index_dict

    def dual_gene(self) -> Dict[int, List[int]]:
        D: Dict[int, List[List[int]]] = {}
        final: Dict[int, List[int]] = {}
        for gr in self.genomic_region_to_SNPs:
            snps = self.genomic_region_to_SNPs[gr]
            start = self.comps
            D[gr] = {s: [] for s in snps} #S
            start = gr
            if len(self.comps[start]) == 1:
                gene = self.genes[list(self.comps[start])[0]]
                final[min(gene.snps)] = sorted(gene.snps)
            else:
                for gene in self.comps[start]:
                    for snp in self.genes[gene].snps:
                        D[gr][snp].append(gene)

                for snp in D[gr]:
                    D[gr][snp] = sorted(D[gr][snp])
                tmp: Dict[List[int], List[int]] = {}
                for s in D[gr]:
                    if D[gr][s] not in tmp:
                        tmp[D[gr][s]] = [s]
                    else:
                        tmp[D[gr][s]].append(s)
                for snps in tmp.values():
                    final[min(snps)] = sorted(snps)
        return final

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

    def make_components(self, snps_to_use: Dict[int, List[int]]):
        comp_mins: List[int] = []
        comps: Dict[int, List[int]] = {}
        for start in snps_to_use:
            snps = snps_to_use[start]
            if len(snps) > 1:
                for snp in snps:
                    comps[snp] = snps
                comp_mins.append(start)
        self.components, self.comp_mins = comps, sorted(comp_mins)

    def find_SNPs(self, g: Gene):
        s = self.find_next_SNP(g, g.bp_start)

        if s == -1:
            g.snps: Set[int] = {}
        else:
            S = self.nodes[s]
            num = S.num
            ALL_SNPS: List[Node] = []
            while S.pos < g.bp_end:
                ALL_SNPS.append(S)
                num += 1
                if num < len(self.nodekeys):
                    S = self.nodes[self.nodekeys[num]]
                else:
                    break
            L = len(ALL_SNPS)
            r_i = 0
            s_i = 0

            snps = {i: [] for i in range(g.exon_count)} #S
            while s_i < L:
                S = ALL_SNPS[s_i]
                r = g.ranges[r_i]
                tf = (S.pos >= r[0], S.pos < r[1])
                if tf[0] and tf[1]:
                    snps[r_i].append(S)
                    s_i += 1
                elif tf[0]:
                    r_i += 1
                elif tf[1]:
                    s_i += 1
            g.snps = {y.index for x in snps.values() for y in x}

    def find_next_SNP(self, g: Gene, b: int) -> int:
        if g.chrom not in self.phasable_positions:
            return -1
        else:
            positions = self.phasable_positions[g.chrom]
            h = bisect_left(positions, b)
            return self.PSdict[g.chrom][positions[h]] if h < len(positions) else -1


def make_data_from_RNA_data(RD: RNAData) -> Data:
    read_list = RD.long_read_list
    data = edges_from_readlist(read_list)
    return Data(
        data, RD.states, RD.k, RD.error, read_list, RD.positions, RD.names, RD.chroms
    )
