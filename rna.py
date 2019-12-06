from graph import Node
from gene import make_genomic_graph, assign_reads_to_genomic_regions
from read import Read
from common import score
from rates import find_rates


def in_range(S: Node, r: tuple[int, int]) -> tuple[bool, bool]:
    p = S.position
    return (p >= r[0], p < r[1])


class RNAData:
    states: dict[int, int]
    initial_genes: dict[int, Gene]
    genes: dict[int, Gene]
    error: float
    all_reads: dict[int, Read]
    positions: dict[int, int]
    names: dict[int, str]
    chroms: dict[int, str]
    isodict: dict[str, list[float, str]]
    k: int

    read_list: dict[int, Read]
    short_read_list: dict[int, Read]
    long_read_list: dict[int, Read]

    nodekeys: list[int]
    nodes: dict[int, Node]
    chrom_list: list[str]
    phasable_positions: dict[str, list[int]]
    PSdict: dict[str, dict[int, int]]

    comps: dict[int, set[int]]
    comps_dict: dict[int, int]
    comps_mins_dict: dict[int, set[int]]
    SNP_to_genomic_region: dict[int, int]
    genomic_region_to_SNPs: dict[int, list[int]]
    reads_by_GR: dict[str, list[Read]]
    reads_by_indiv_gene: dict[int, dict[int, Read]]
    final: dict[int, list[int]]
    snps_to_use: dict[int, list[int]]

    reads_by_comsnps: dict[int, list[Read]]
    read_dict: dict[Read, list[Read]]
    counts: dict[Read, list[int]]
    LLsnp: dict[int, float]
    rates: dict[int, dict[int, float]]

    components: dict[int, list[int]]
    comp_mins: list[int]


    def __init__(
        self: RNAData,
        states: dict[int, int],
        genes: dict[int, Gene],
        filtered_genes: dict[int, Gene],
        error: float,
        read_list: dict[int, Read],
        positions: dict[int, int],
        names: dict[int, str],
        chroms: dict[int, str],
        isodict: dict[str, list[float, str]]
    ):
        self.k = 2
        self.states = states
        self.initial_genes = genes
        self.genes = filtered_genes
        self.error = error
        self.all_reads = read_list
        self.positions = positions
        self.names = names
        self.chroms = chroms
        self.isodict = isodict

        short_read_list = dict[int, Read]()
        long_read_list = dict[int, Read]()
        for i in self.all_reads:
            r = self.all_reads[i]
            r.rates = {0: 0.5, 1: 0.5}
            if len(r.keys) == 1:
                short_read_list[i] = r
            else:
                long_read_list[i] = r

        self.read_list = short_read_list
        self.short_read_list = short_read_list
        self.long_read_list = long_read_list

        self.nodekeys = sorted(list(self.node_keys()))

        d = dict[int, Node]()
        chrom_list = dict[str, int]()
        for num in range(len(self.nodekeys)):
            index = self.nodekeys[num]
            d[index] = Node(
                self.k,
                index,
                0,
                [],
                self.positions[index],
                self.names[index],
                self.chroms[index],
                num,
            )
            chrom_list[self.chroms[index]] = 0
        self.nodes = d
        self.chrom_list = sorted(chrom_list.keys())
        self.phasable_positions = {chrom: list[int]() for chrom in self.chrom_list}
        for x in self.nodekeys:
            self.phasable_positions[self.chroms[x]].append(self.positions[x])
        self.PSdict = self.make_position_dict()

        ######################################################################

        print("running through genes...")
        for g in self.genes:
            self.find_SNPs(self.genes[g])
        make_genomic_graph(self.genes)

        ######################################################################

        reads = self.all_reads
        (
            self.comps,
            self.comps_dict,
            self.comps_mins_dict,
            self.SNP_to_genomic_region,
            self.genomic_region_to_SNPs,
            self.reads_by_GR,
        ) = assign_reads_to_genomic_regions(self.genes, reads)
        self.reads_by_indiv_gene = self.assign_reads_to_indiv_genes()

        ######################################################################

        # Here is where we decide which types of genes to use, those with no isoforms ("no_splicing")
        # all SNPs such that if there are multiple isoforms, those SNPs fall into all of them ("final")

        self.final = self.dual_gene()
        self.snps_to_use = self.find_snps_to_use_all()

        ######################################################################

        self.reads_by_comsnps = self.find_reads_by_comsnps()
        self.read_dict = self.make_read_dict()  ##1-reads only
        self.counts = self.assign_counts_to_indiv_snps()
        self.LLsnp = self.assign_LL_dif_to_snps()
        print("assigning rates")
        self.rates = self.assign_rates2(self.snps_to_use)

    ######################################################################

    def node_keys(self: RNAData) -> list[int]:
        s = {k: 0 for r in self.read_list.values() for k in r.keys}
        return list(s.keys())

    def make_position_dict(self: RNAData) -> dict[str, dict[int, int]]:
        position_to_SNP_index_dict = {chrom: dict[int, int]() for chrom in self.chrom_list}
        for s in self.nodekeys:
            S = self.nodes[s]
            c = S.chrom
            p = S.position
            position_to_SNP_index_dict[c][p] = S.index
        return position_to_SNP_index_dict

    def assign_reads_to_indiv_genes(self: RNAData) -> dict[int, dict[int, Read]]:
        D = {gr: dict[int, Read]() for gr in self.comps_mins_dict}
        for gr in self.comps_mins_dict:
            for g in self.comps[gr]:
                D[gr][g] = list[]()
                snps = self.genes[g].snps
                for r in self.reads_by_GR[gr]:
                    tf = True
                    for k in r.keys:
                        tf = tf and (k in snps)
                    if tf:
                        D[gr][g].append(r)
        return D

    def dual_gene(self: RNAData) -> dict[int, list[int]]:
        D = dict[int, list[list[int]]]()
        final = dict[int, list[int]]()
        for gr in self.genomic_region_to_SNPs:
            snps = self.genomic_region_to_SNPs[gr]
            start = self.comps
            D[gr] = {s: list[list[int]]() for s in snps}
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
                tmp = dict[list[int], list[int]]()
                for s in D[gr]:
                    if D[gr][s] not in tmp:
                        tmp[D[gr][s]] = [s]
                    else:
                        tmp[D[gr][s]].append(s)
                for snps in tmp.values():
                    final[min(snps)] = sorted(snps)
        return final

    def find_snps_to_use_all(self: RNAData) -> dict[int, list[int]]:
        snps_to_use = dict[int, list[int]]()
        for start in self.genomic_region_to_SNPs:
            snps = self.genomic_region_to_SNPs[start]
            if len(snps) > 1:
                snps_to_use[min(snps)] = sorted(snps)
        return snps_to_use

    def find_reads_by_comsnps(self: RNAData) -> dict[int, list[Read]]:
        # at this point we will only use reads that fall strictly within common snps
        # we should see how many long reads we arent using.
        # we should consider adding those back in to make components (earlier)
        reads_by_comsnps = dict[int, list[Read]]()
        for start in sorted(self.snps_to_use.keys()):
            comp = self.snps_to_use[start]
            m = start
            reads_by_comsnps[start] = list[Read]()
            gr = self.SNP_to_genomic_region[start]
            for r in self.reads_by_GR[gr]:
                tf = True
                for k in r.keys:
                    tf = tf and (k in comp)
                if tf:
                    reads_by_comsnps[start].append(r)
        return reads_by_comsnps

    def make_read_dict(self: RNAData) -> dict[Read, list[Read]]:
        read_dicts = dict[Read, list[Read]]()
        for m in self.reads_by_comsnps:
            for r in self.reads_by_comsnps[m]:
                if len(r.keys) == 1:
                    k = r.keys[0]
                    if k in read_dicts:
                        read_dicts[k].append(r)
                    else:
                        read_dicts[k] = [r]
                else:
                    for k in r.keys:
                        new_read = Read({k: r.read[k]}, r.count, -1)
                        if k in read_dicts:
                            read_dicts[k].append(new_read)
                        else:
                            read_dicts[k] = [new_read]
        return read_dicts

    def assign_counts_to_indiv_snps(self: RNAData) -> dict[Read, list[int]]:
        ###ASSUMES 1READS ONLY #1reads
        all_counts = dict[Read, list[int]]()
        for s in self.read_dict:
            counts = [0, 0]
            for r in self.read_dict[s]:
                counts[r.read[s] % 2] += r.count
            all_counts[s] = counts
        return all_counts

    def assign_LL_dif_to_snps(self: RNAData) -> dict[int, float]():
        LL = dict[int, float]()
        for start in self.snps_to_use:
            for snp in self.snps_to_use[start]:
                LL[snp] = score(self.counts[snp])
        return LL  # LL_gr,LL

    def assign_rates2(
        self: RNAData, snps_to_use: dict[int, list[int]]
    ) -> dict[int, dict[int, float]]:
        #times = {}
        rates = dict[int, dict[int, float]]()
        for start in snps_to_use:
            snps = sorted(snps_to_use[start])
            reads = [r for s in snps for r in self.read_dict[s]]
            # t = time.time()
            r = find_rates(snps, reads, 0.6)  ##Choose adjacent snp rate
            # t = t - time.time()
            # times[t] = snps
            for s in snps:
                rates[s] = r
            for read in reads:
                read.rates = r
        #self.times = times
        return rates

    ############################################################

    def make_components(self: RNAData, snps_to_use: dict[int, list[int]]):
        comp_mins = list[int]()
        comps = dict[int, list[int]]()
        for start in snps_to_use:
            snps = snps_to_use[start]
            if len(snps) > 1:
                for snp in snps:
                    comps[snp] = snps
                comp_mins.append(start)
        self.components, self.comp_mins = comps, sorted(comp_mins)

    def find_SNPs(self: RNAData, g: Gene):
        s = self.find_next_SNP(g, g.bp_start)

        if s == -1:
            g.snps = set[int]()
        else:
            S = self.nodes[s]
            num = S.num
            ALL_SNPS = list[Node]()
            while S.position < g.bp_end:
                ALL_SNPS.append(S)
                num += 1
                if num < len(self.nodekeys):
                    S = self.nodes[self.nodekeys[num]]
                else:
                    break
            L = len(ALL_SNPS)
            r_i = 0
            s_i = 0

            snps = {i: [] for i in range(g.exon_count)}

            while s_i < L:
                S = ALL_SNPS[s_i]
                r = g.ranges[r_i]
                tf = in_range(S, r)
                if tf[0] and tf[1]:
                    snps[r_i].append(S)
                    s_i += 1
                elif tf[0]:
                    r_i += 1
                elif tf[1]:
                    s_i += 1
            g.snps = {y.index for x in snps.values() for y in x}

    def find_next_SNP(self: RNAData, g: Gene, b: int) -> int:
        if g.chrom not in self.phasable_positions:
            return -1
        else:
            positions = self.phasable_positions[g.chrom]
            L = len(positions)
            PSdict = self.PSdict[g.chrom]
            nodekeys = self.nodekeys
            h = bisect.bisect_left(positions, b)
            if h < L:
                i = positions[h]
                j = PSdict[i]
                return j
            else:
                return -1


def make_data_from_RNA_data(RD: RNAData) -> Data:
    read_list = RD.long_read_list
    data = edges_from_readlist(read_list)
    return DATA(
        data, RD.states, RD.k, RD.error, read_list, RD.positions, RD.names, RD.chroms
    )