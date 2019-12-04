import alg
import RD_updates
import global_vars

##import prescoring
import time
import copy

from monkeytype.encoding import DUMMY_NAME
from rna_class import RNA_DATA
from typing import Any, Dict, List, Tuple


def stats(
    RD: RNA_DATA,
    size_factor: int,
    rate_factor: float,
    rate_cutoff: float,
    coverage_cutoff: int,
    rate_dep_cutoff: int,
    cutoff: float,
    conf: float,
    show: bool = False,
) -> None:

    # print "RD filtering parameters"
    # print "size_factor: ", size_factor
    # print "rate_factor: ", rate_factor
    # print "rate_cutoff: ", rate_cutoff
    # print "coverage_cutoff: ", coverage_cutoff
    # print "rate_dep_cutoff: ", rate_dep_cutoff
    # print "cutoff: ", cutoff
    # print "conf: ", conf

    ### all of these current_STU* filter the set of SNPs we are going to try to phase
    current_STU0 = RD.snps_to_use
    print("Original RD SNP dictionary size: ", len(current_STU0), num_SNP(current_STU0))
    current_STU1 = RD_updates.update_snps_to_use_rate(RD, current_STU0, rate_cutoff)
    print(
        "STU1(rate_cutoff) RD SNP dictionary size: ",
        len(current_STU1),
        num_SNP(current_STU1),
    )
    current_STU2 = RD_updates.update_snps_to_use_coverage(
        RD, current_STU1, coverage_cutoff
    )
    print(
        "STU2(coverage_cutoff) RD SNP dictionary size: ",
        len(current_STU2),
        num_SNP(current_STU2),
    )
    current_STU3 = RD_updates.update_snps_to_use_size_cluster(
        RD, current_STU2, size_factor
    )
    print(
        "STU3(size_factor) RD SNP dictionary size: ",
        len(current_STU3),
        num_SNP(current_STU3),
    )
    current_STU4 = RD_updates.update_snps_to_use_rate_dependent_cutoff(
        RD, current_STU3, rate_dep_cutoff, conf
    )
    print(
        "STU4(rate_dep_cutoff,conf) RD SNP dictionary size: ",
        len(current_STU4),
        num_SNP(current_STU4),
    )
    #    print current_STU4
    current_STU5 = RD_updates.update_snps_to_use_cutoff(RD, current_STU4, cutoff)
    print(
        "STU5(cutoff) RD SNP dictionary size: ",
        len(current_STU5),
        num_SNP(current_STU5),
    )
    #    print current_STU5
    current_STU6 = RD_updates.update_snps_to_use_rate_cluster(
        RD, current_STU5, rate_factor
    )
    print(
        "STU6(rate_factor) RD SNP dictionary size: ",
        len(current_STU6),
        num_SNP(current_STU6),
    )
    #    print current_STU6

    RD.make_components(current_STU6)
    X = alg.RNA_phase(
        0.0001,
        global_vars.pair_thresh,
        RD.error,
        RD.read_dict,
        RD.comp_mins,
        RD.components,
    )
    S = [switches_in_comp2(x, X, global_vars.V, RD)[0] for x in X]
    # print sum(S), len(RD.components), len(X)
    if len(RD.components) == 0:
        print(
            "Warning: No fragments covering only 1 SNPs were found in the RNAfragmat input"
        )
        print(
            "HapTree-X will not consider differential allele-specific expression (DASE) during phasing"
        )
        print(
            "If you would like to use DASE as well, please run chair with parameter 1"
        )

    C = con_dis_non(X, global_vars.V)[1]
    # print sum(C),C[0],C[1]
    # if sum(C) != 0:
    #    print C[1]/float(sum(C))
    # if C[0]+C[1] != 0:
    #    print C[1]/float(C[0]+C[1])
    if show:
        for x in X:
            if switches_in_comp2(x, X, global_vars.V, RD)[0] > 0:
                gr = RD.SNP_to_genomic_region[x]
                # print RD.genes[gr].gene_type
                # print x, RD.components[x]
                # print get_counts2(RD.components[x],RD.read_dict)
                # for y in sorted(X[x][0]):
                #                print X[x][0][y],
                # for y in sorted(X[x][0]):
                #                print global_vars.V[0][y],


def num_SNP(STU: Dict[int, List[int]]) -> int:
    return sum(map(len, list(STU.values())))


def get_counts2(snps, read_dict):
    ##requires all reads to have length 1
    reads = []
    for s in snps:
        for r in read_dict[s]:
            reads.append(r)

    counts = {i: [0, 0] for i in range(len(snps))}
    back = {}
    for i in range(len(snps)):
        back[snps[i]] = i
    for R in reads:
        for snp in R.keys:
            counts[back[snp]][R.read[snp]] += R.count
    return counts


def comp_phase_print(RD, snps_to_use):
    RD.make_components(snps_to_use)
    X = alg.RNA_phase(RD, 0.0001, global_vars.pair_thresh)
    S = [switches_in_comp2(x, X, global_vars.V, RD)[0] for x in X]
    # print sum(S), len(RD.components), sum(S)/float(len(RD.components)),len(RD.components),


def con_dis_non(
    X: Dict[int, Dict[int, Dict[int, int]]], V: Dict[int, Dict[int, str]]
) -> Tuple[Dict[int, DUMMY_NAME], List[int]]:
    counts = {}  # {x:{} for x in X}
    for x in X:
        countsx = {0: 0, "n": 0, 1: 0}
        for y in X[x][0]:
            if V[0][y] == ".":
                countsx["n"] += 1
            else:
                if V[0][y] == X[x][0][y]:
                    countsx[0] += 1
                else:
                    countsx[1] += 1
        m, M = sorted([countsx[0], countsx[1]])
        n = countsx["n"]
        counts[x] = {"c": M, "n": n, "d": m}
    totals_con = sum([counts[x]["c"] for x in X])
    totals_dis = sum([counts[x]["d"] for x in X])
    totals_none = sum([counts[x]["n"] for x in X])
    totals_all = [totals_con, totals_dis, totals_none]
    total = sum(totals_all)
    return counts, totals_all


def m_swcounter2(sol, gold, G):
    switches = 0
    for start in G.comp_mins:
        num, l = switches_in_comp2(start, sol, gold, G)
        switches += num

    return switches


def switches_in_comp2(
    start: int,
    sol: Dict[int, Dict[int, Dict[int, int]]],
    gold: Dict[int, Dict[int, str]],
    G: RNA_DATA,
) -> Tuple[int, List[Any]]:
    m0, switches0 = switches_comp_strand2(start, 0, sol, gold, G)
    m1, switches1 = switches_comp_strand2(start, 1, sol, gold, G)
    if m0 < m1:
        return m0, switches0
    else:
        return m1, switches1


def switches_comp_strand2(
    start: int,
    strand: int,
    sol: Dict[int, Dict[int, Dict[int, int]]],
    gold: Dict[int, Dict[int, str]],
    G: RNA_DATA,
) -> Tuple[int, List[Any]]:
    m = 0
    switches = []
    for i in G.components[start]:
        if sol[start][strand][i] == gold[0][i] or gold[0][i] == ".":
            None
        else:
            switches.append(i)
            m += 1
            strand = 1 - strand
    return m, switches


def switches_comp_strand3(start, strand, sol, gold):
    m = 0
    switches = []
    for i in sol[start][0].keys():
        if sol[start][strand][i] == gold[0][i] or gold[0][i] == ".":
            None
        else:
            switches.append(i)
            m += 1
            strand = 1 - strand
    return m, switches


def switches_in_comp3(start, sol, gold):
    m0, switches0 = switches_comp_strand3(start, 0, sol, gold)
    m1, switches1 = switches_comp_strand3(start, 1, sol, gold)
    if m0 < m1:
        return m0, switches0
    else:
        return m1, switches1


def restrict_phasing(snp_to_gr, phasing_comps_old):
    # gr is "genomic_region", cluster of overlapping genes
    phasing_comps_new = {}
    for m in phasing_comps_old:
        local_comps_by_gr = {}
        comp = phasing_comps_old[m][0]
        for snp in comp:
            if snp in snp_to_gr:
                gr = snp_to_gr[snp]
                if gr in local_comps_by_gr:
                    local_comps_by_gr[gr].append(snp)
                else:
                    local_comps_by_gr[gr] = [snp]
        for gr in local_comps_by_gr:
            sub_snps_gr = local_comps_by_gr[gr]
            if len(sub_snps_gr) > 1:
                comp_min = min(local_comps_by_gr[gr])
                comp_phase = {0: {}, 1: {}}
                for snp in sub_snps_gr:
                    comp_phase[0][snp] = phasing_comps_old[m][0][snp]
                    comp_phase[1][snp] = phasing_comps_old[m][1][snp]
                phasing_comps_new[comp_min] = comp_phase

    return phasing_comps_new


def phasing_completeness(phasing):
    num_comps = len(phasing)
    num_snps = sum([len(x[0]) for x in phasing.values()])
    print(
        "num_comps: ",
        num_comps,
        "num_snps: ",
        num_snps,
        "completeness:  ",
        num_snps - num_comps,
        end=" ",
    )
    return num_snps - num_comps


def count_all_switches(phasing, V):
    S = [switches_in_comp3(x, phasing, V)[0] for x in phasing]
    return sum(S)
