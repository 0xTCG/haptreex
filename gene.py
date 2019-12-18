from graph import Node
from read import Read
from typing import Tuple, Dict, List, Set
from dataclasses import dataclass


NOT_SEEN = -1
OVER_SEEN = -2



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
    genes: List[Gene], reads: List[Read]
):
    comps, comps_dict = build_all_genomic_comps(genes)
    pprint({ genes[g].name: sorted(genes[q].name for q in y) for g, y in self.comps.items()})
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
    for j in genes:
        g = genes[j]
        tid = g.name
        if tid in isodict:
            val = isodict[tid][0]
            if val > 0:
                if val / gene_cov[isodict[tid][1]] > 0.01:
                    new_genes[j] = g
    return new_genes

