from bisect import bisect_left
from graph import build_components, Component
from read import Read
from rates import find_rates, score
from typing import Tuple, Dict, List, Set
from dataclasses import dataclass
from math import log


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
        self.exons = exons

        self.neighbors: Set[int] = set()
        self.snps: Set[int] = set()

    def __eq__(self, other):
        return self.name == other.name

    def set_snps(self, positions: List[Tuple[int, int]]):
        h = bisect_left(positions, (self.interval[0], 0))
        if h == len(positions):
            return
        if not self.exons:
            return
        exon = 0
        self.snps = set()
        for p in range(h, len(positions)):
            pos, snp = positions[p]
            pos += 1
            while exon < len(self.exons) and pos >= self.exons[exon][1]:
                exon += 1
            if pos >= self.interval[1] or exon >= len(self.exons):
                break
            if pos >= self.exons[exon][0]:
                if pos < self.exons[exon][1]:
                    self.snps.add(snp)
                else:
                    exon += 1


@dataclass
class RNAGraph:
    genes: List[Gene]
    ploidy: int
    multi_reads: List[Read]
    counts: Dict[int, List[int]]
    LL: Dict[int, float]

    components: Dict[int, Component]  # {Component ID: Component}
    component_index: Dict[int, int]  # SNP ID: Component ID
    snp_reads: Dict[int, List[Read]]  # Reads for each SNP

    def __init__(
        self,
        vcf,
        genes: List[Gene],
        reads: List[Read],
        size_factor: int = 2,
        rate_factor: float = 0.2,
        rate_cutoff: float = 0.6,
        coverage_cutoff: int = 0,
        rate_dep_cutoff: int = 2,
        cutoff: float = 0.001,
        conf: float = 0.2
    ):
        self.ploidy = 2
        self.multi_reads = [r for r in reads if len(r.snps) > 1]

        # Initialize gene nodes
        # (1) Initialize gene SNPs
        phasable_positions: Dict[str, List[Tuple[int, int]]] = {}
        # Add `if len(r.snps) >= 1` to the set below to enable the original behaviour
        for s in {k for r in reads for k in r.snps}:
            snp = vcf.snps[s]
            phasable_positions.setdefault(snp.chr, []).append((snp.pos, snp.id))
        for p in phasable_positions.values():
            p.sort()
        for g in genes:
            if g.chr in phasable_positions:
                g.set_snps(phasable_positions[g.chr])
        # (2) Initialize gene neighbours (genes that share a common SNP)
        snp_to_genes: Dict[int, Set[int]] = {}  # SNP to gene
        for gi, g in enumerate(genes):
            for s in g.snps:
                snp_to_genes.setdefault(s, set()).add(gi)
        for g in genes:
            for snp in g.snps:
                for n in snp_to_genes[snp]:
                    g.neighbors.add(n)

        # Build gene components from genes that have at least one SNP
        gene_components, _ = build_components(
            dict((i, g) for i, g in enumerate(genes) if g.snps)
        )
        snp_to_comp = {  # SNP -> genic component
            snp: root
            for root, comp in gene_components.items()
            for gene in comp.nodes
            for snp in genes[gene].snps
        }
        # self.final = self.dual_gene()  # Which types of genes to use?
        comp_to_snps: Dict[int, List[int]] = {m: [] for m in gene_components}
        for s in snp_to_comp:
            comp_to_snps[snp_to_comp[s]].append(s)
        comp_snps: Dict[int, List[int]] = {}  # root SNP -> list of SNPs in component
        for root, snps in comp_to_snps.items():
            if len(snps) > 1:
                snps.sort()
                comp_snps[snps[0]] = sorted(snps)
        # Only use reads that fall strictly within a gene
        comp_to_reads: Dict[int, List[Read]] = {root: [] for root in gene_components}
        NOTSEEN = -1  # NOTSEEN span non-exonic regions
        OVERSEEN = -2  # OVERSEEN typically span exonic and intronic regions
        comp_to_reads[OVERSEEN], comp_to_reads[NOTSEEN] = [], []
        for read in reads:
            comps = {snp_to_comp.get(snp, NOTSEEN) for snp in read.snps}
            # comps -= {NOTSEEN}  # Uncomment to allow reads with partial exonic span
            comp = comps.pop() if len(comps) == 1 else OVERSEEN
            assert comp in comp_to_reads
            comp_to_reads[comp].append(read)
        self.snp_reads: Dict[int, List[Read]] = {}  # SNP -> 1-reads in a component
        for root in comp_snps:
            for read in comp_to_reads[snp_to_comp[root]]:
                if len(read.snps) == 1:
                    snp = list(read.snps.keys())[0]
                    self.snp_reads.setdefault(snp, []).append(read)
                else:
                    for snp in read.snps:
                        self.snp_reads.setdefault(snp, []).append(
                            Read({snp: read.snps[snp]}, read.count, -1)
                        )
        for root, snps in comp_snps.items():
            # Ignore SNPs that are not expressed by reads that survived till here
            comp_snps[root] = [s for s in snps if s in self.snp_reads]
            # print(f'SNP {snp} not expressed in reads!')
        
        # Calculate counts and LL_dif for filtering
        self.counts: Dict[int, List[int]] = {}
        for snp in self.snp_reads:
            self.counts[snp] = [0, 0]
            for read in self.snp_reads[snp]:
                self.counts[snp][read.snps[snp] % 2] += read.count
        self.LL = {s: score(self.counts[s]) for ss in comp_snps.values() for s in ss}

        # Assign DASE rates
        self.rates: Dict[int, List[float]] = {}
        for snps in comp_snps.values():
            reads = [read for snp in snps for read in self.snp_reads[snp]]
            rate = find_rates(snps, reads, 0.6)  # Choose adjacent snp rate
            for snp in snps:
                self.rates[snp] = rate
            for read in reads:
                read.rates = rate

        def num(ss):
            return sum(len(s) for _, s in ss.items())

        # Apply filters
        print("Original RD SNP dictionary size: ", len(comp_snps), num(comp_snps))

        comp_snps = self.rate_filter(comp_snps, rate_cutoff)
        print("STU1(rate_cutoff) RD SNP dictionary size: ", len(comp_snps), num(comp_snps))

        comp_snps = self.coverage_filter(comp_snps, coverage_cutoff)
        print("STU2(coverage_cutoff) RD SNP dictionary size: ", len(comp_snps), num(comp_snps))

        comp_snps = self.size_cluster_filter(comp_snps, size_factor)
        print("STU3(size_factor) RD SNP dictionary size: ", len(comp_snps), num(comp_snps))

        comp_snps = self.rate_dependent_filter(comp_snps, rate_dep_cutoff, conf)
        print("STU4(rate_dep_cutoff,conf) RD SNP dictionary size: ", len(comp_snps), num(comp_snps))

        comp_snps = self.cutoff_filter(comp_snps, cutoff)
        print("STU5(cutoff) RD SNP dictionary size: ", len(comp_snps), num(comp_snps))

        comp_snps = self.rate_cluster_filter(comp_snps, rate_factor)
        print("STU6(rate_factor) RD SNP dictionary size: ", len(comp_snps), num(comp_snps))

        # Make SNP connected components
        self.components = {}
        self.component_index = {}
        for root in comp_snps:
            snps = comp_snps[root]
            if len(snps) > 1:
                self.components[root] = Component(root, snps, [])
                for snp in snps:
                    self.component_index[snp] = root

    def rate_filter(self, snps: Dict[int, List[int]], rate_cutoff: float):
        # High-confidence SNPs
        return {r: s for r, s in snps.items() if max(self.rates[r]) >= rate_cutoff}

    def coverage_filter(self, snps: Dict[int, List[int]], coverage_cutoff: int):
        new: Dict[int, List[int]] = {}
        for r in snps:
            covered_snps = sorted([
                snp for snp in snps[r] if sum(self.counts[snp]) > coverage_cutoff
            ])
            if len(covered_snps) > 1:
                new[covered_snps[0]] = covered_snps
        return new

    def size_cluster_filter(self, snps: Dict[int, List[int]], size_factor: int):
        return {
            min(cluster): sorted(cluster)
            for r in snps
            for cluster in self.clusters_size(snps[r], size_factor)
        }

    def rate_dependent_filter(
        self, snps: Dict[int, List[int]], cutoff: int, conf: float
    ):
        new: Dict[int, List[int]] = {}
        for r in snps:
            confident_snps = sorted([
                snp for snp in snps[r]
                if self.new_score(snp, conf) > cutoff
            ])
            if len(confident_snps) > 1:
                new[confident_snps[0]] = confident_snps
        return new

    def cutoff_filter(self, snps: Dict[int, List[int]], cutoff: float):
        new: Dict[int, List[int]] = {}
        for r in snps:
            confident_snps = sorted([
                snp for snp in snps[r] if self.LL[snp] < cutoff
            ])
            if len(confident_snps) > 1:
                new[confident_snps[0]] = confident_snps
        return new

    def rate_cluster_filter(self, snps: Dict[int, List[int]], rate_factor: float):
        return {
            min(cluster): sorted(cluster)
            for r in snps
            for cluster in self.clusters_rate(snps[r], rate_factor)
            if len(cluster) > 0
        }

    def clusters_size(self, snps: List[int], sum_factor: int):
        def sort_key(a: Tuple[int, float]) -> float:
            return a[1]
        tmp = {s: float(sum(self.counts[s])) for s in snps}
        order = sorted(tmp.items(), key=sort_key)

        clusters = [[order[0][0]]]
        for j in range(len(order) - 1):
            s0 = order[j][0]
            s1 = order[j + 1][0]
            if tmp[s1] < sum_factor * tmp[s0]:
                clusters[-1].append(s1)
            else:
                clusters.append([s1])
        return clusters

    def clusters_rate(self, snps: List[int], rate_factor: float):
        def sort_key(a: Tuple[int, float]) -> float:
            return a[1]
        tmp = {s: max(self.counts[s]) / float(sum(self.counts[s])) - 0.5 for s in snps}
        order = sorted(tmp.items(), key=sort_key)

        clusters = [[order[0][0]]]
        for j in range(len(order) - 1):
            s0 = order[j][0]
            s1 = order[j + 1][0]
            if tmp[s1] < rate_factor + tmp[s0]:
                clusters[-1].append(s1)
            else:
                clusters.append([s1])
        return clusters

    def new_score(self, snp: int, conf: float):
        pair = self.counts[snp]
        rates = self.rates[snp]
        r = conf * max(min(rates), 0.05) + 0.5 * (1 - conf)
        k, n = min(pair), sum(pair)
        return max(-(n - 2 * k) * log(r / (1 - r)), log(2.0))
