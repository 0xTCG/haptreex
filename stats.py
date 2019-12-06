from alg import RNA_phase
from stats_helper import *
from rna import RNAData
import common


def num_SNP(STU: dict[int, list[int]]) -> int:
    return sum(len(s) for s in list(STU.values()))


def switches_comp_strand2(
    start: int,
    strand: int,
    sol: dict[int, dict[int, dict[int, int]]],
    gold: dict[int, dict[int, str]],
    G: RNAData,
) -> tuple[int, list[int]]:
    m = 0
    switches = list[int]()
    for i in G.components[start]:
        if sol[start][strand][i] != gold[0][i] and gold[0][i] != ".":
            switches.append(i)
            m += 1
            strand = 1 - strand
    return m, switches


def switches_in_comp2(
    start: int,
    sol: dict[int, dict[int, dict[int, int]]],
    gold: dict[int, dict[int, str]],
    G: RNAData,
) -> tuple[int, list[int]]:
    m0, switches0 = switches_comp_strand2(start, 0, sol, gold, G)
    m1, switches1 = switches_comp_strand2(start, 1, sol, gold, G)
    if m0 < m1:
        return m0, switches0
    else:
        return m1, switches1


def con_dis_non(
    X: dict[int, dict[int, dict[int, int]]], V: dict[int, dict[int, str]]
) -> tuple[dict[int, dict[str, int]], list[int]]:
    counts = dict[int, dict[str, int]]()  # {x:{} for x in X}
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


def stats(
    RD: RNAData,
    size_factor: int,
    rate_factor: float,
    rate_cutoff: float,
    coverage_cutoff: int,
    rate_dep_cutoff: int,
    cutoff: float,
    conf: float,
    show: bool = False,
):
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
        common.pair_thresh,
        RD.error,
        RD.read_dict,
        RD.comp_mins,
        RD.components,
    )
    S = [switches_in_comp2(x, X, common.V, RD)[0] for x in X]
    # print sum(S), len(RD.components), len(X)
    if len(RD.components) == 0:
        print(
            "Warning: No fragments covering only 1 SNPs were found in the RNAfragmat "
            + "input.\nHapTree-X will not consider differential allele-specific "
            + "expression (DASE) during phasing.\nIf you would like to use DASE as "
            + "well, please run chair with parameter 1"
        )

    C = con_dis_non(X, common.V)[1]
    if show:
        for x in X:
            if switches_in_comp2(x, X, common.V, RD)[0] > 0:
                gr = RD.SNP_to_genomic_region[x]


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
