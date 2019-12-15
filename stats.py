from alg import RNA_phase
from stats_helper import *
from rna import RNAData
import common
from typing import Tuple, Dict, List, Set


def num_SNP(STU: Dict[int, List[int]]) -> int:
    return sum(len(s) for s in STU.values())


def switches_comp_strand2(
    start: int,
    strand: int,
    sol: Dict[int, Tuple[Dict[int, int], Dict[int, int]]],
    gold: Dict[int, Dict[int, int]],
    G: RNAData
) -> Tuple[int, List[int]]:
    m = 0
    switches: List[int] = []
    for i in G.components[start]:
        p = [sol[start][0][i], sol[start][1][i]]
        if p[strand] != gold[0][i] and gold[0][i] != common.DOT:
            switches.append(i)
            m += 1
            strand = 1 - strand
    return m, switches


def switches_in_comp2(
    start: int,
    sol: Dict[int, Tuple[Dict[int, int], Dict[int, int]]],
    gold: Dict[int, Dict[int, int]],
    G: RNAData
) -> Tuple[int, List[int]]:
    m0, switches0 = switches_comp_strand2(start, 0, sol, gold, G)
    m1, switches1 = switches_comp_strand2(start, 1, sol, gold, G)
    if m0 < m1:
        return m0, switches0
    else:
        return m1, switches1


def con_dis_non(
    X: Dict[int, Tuple[Dict[int, int], Dict[int, int]]], V: Dict[int, Dict[int, int]]
) -> Tuple[Dict[int, Dict[str, int]], List[int]]:
    counts: Dict[int, Dict[str, int]] = {}  # {x:{} for x in X}
    for x in X:
        countsx = {0: 0, -1: 0, 1: 0}
        for y in X[x][0]:
            if V[0][y] == common.DOT:
                countsx[-1] += 1
            else:
                if V[0][y] == X[x][0][y]:
                    countsx[0] += 1
                else:
                    countsx[1] += 1
        m, M = sorted([countsx[0], countsx[1]])
        n = countsx[-1]
        counts[x] = {"c": M, "n": n, "d": m}
    totals_con = sum([counts[x]["c"] for x in X])
    totals_dis = sum([counts[x]["d"] for x in X])
    totals_none = sum([counts[x]["n"] for x in X])
    totals_all = [totals_con, totals_dis, totals_none]
    total = sum(totals_all)
    return counts, totals_all


def stats(
    RD: RNAData,
    size_factor: int,
    rate_factor: float,
    rate_cutoff: float,
    coverage_cutoff: int,
    rate_dep_cutoff: int,
    cutoff: float,
    conf: float,
    show: bool = False
):
    ### all of these current_STU* filter the Set of SNPs we are going to try to phase
    current_STU0 = RD.snps_to_use
    print(f"Original RD SNP dictionary size: {len(current_STU0)} {num_SNP(current_STU0)}")
    #for s in sorted(current_STU0):
    #    print(f"{s} {max(RD.rates[s].values())} {rate_cutoff}")
    current_STU1 = update_snps_to_use_rate(RD, current_STU0, rate_cutoff)
    print(f"STU1(rate_cutoff) RD SNP size: {len(current_STU1)} {num_SNP(current_STU1)}")
    current_STU2 = update_snps_to_use_coverage(RD, current_STU1, coverage_cutoff)
    print(f"STU2(coverage_cutoff) RD SNP size: {len(current_STU2)} {num_SNP(current_STU2)}")
    current_STU3 = update_snps_to_use_size_cluster(RD, current_STU2, size_factor)
    print(f"STU3(size_factor) RD SNP size: {len(current_STU3)} {num_SNP(current_STU3)}")
    current_STU4 = update_snps_to_use_rate_dependent_cutoff(RD, current_STU3, rate_dep_cutoff, conf)
    print(f"STU4(rate_dep_cutoff,conf) RD SNP size: {len(current_STU4)} {num_SNP(current_STU4)}")
    current_STU5 = update_snps_to_use_cutoff(RD, current_STU4, cutoff)
    print(f"STU5(cutoff) RD SNP size: {len(current_STU5)} {num_SNP(current_STU5)}")
    current_STU6 = update_snps_to_use_rate_cluster(RD, current_STU5, rate_factor)
    print(f"STU6(rate_factor) RD SNP size: {len(current_STU6)} {num_SNP(current_STU6)}")
    
    RD.make_components(current_STU6)
    X = RNA_phase(
        0.0001,
        common.PAIR_TRESH,
        RD.error,
        RD.read_dict,
        RD.comp_mins,
        RD.components
    )
    S = [switches_in_comp2(x, X, common.V, RD)[0] for x in X]
    if len(RD.components) == 0:
        print (
            f"Warning: No fragments covering only 1 SNPs were found in the RNAfragmat "
            + f"input.\nHapTree-X will not consider differential allele-specific "
            + f"expression (DASE) during phasing.\nIf you would like to use DASE as "
            + f"well, please run chair with parameter 1"
        )

    C = con_dis_non(X, common.V)[1]
    if show:
        for x in X:
            if switches_in_comp2(x, X, common.V, RD)[0] > 0:
                gr = RD.SNP_to_genomic_region[x]
