from math import log
from rna import RNAData
from typing import Tuple, Dict, List, Set


def sort_key(a: Tuple[int, float]) -> float:
    return a[1]


def clusters_size(
    snps: List[int], counts: Dict[int, List[int]], sum_factor: int
) -> List[List[int]]:
    tmp_dict = {s: float(sum(counts[s])) for s in snps}
    order = sorted(tmp_dict.items(), key=sort_key)
    
    clusters = [[order[0][0]]]
    for j in range(len(order) - 1):
        s0 = order[j][0]
        s1 = order[j + 1][0]
        if tmp_dict[s1] < sum_factor * tmp_dict[s0]:
            clusters[-1].append(s1)
        else:
            clusters.append([s1])
    return clusters


def clusters_rate(
    snps: List[int], counts: Dict[int, List[int]], rate_factor: float
) -> List[List[int]]:
    tmp_dict = {s: max(counts[s]) / float(sum(counts[s])) - 0.5 for s in snps}
    order = sorted(tmp_dict.items(), key=sort_key)

    clusters = [[order[0][0]]]
    for j in range(len(order) - 1):
        s0 = order[j][0]
        s1 = order[j + 1][0]
        if tmp_dict[s1] < rate_factor + tmp_dict[s0]:
            clusters[-1].append(s1)
        else:
            clusters.append([s1])
    return clusters


def new_score(pair: List[int], rates: Tuple[float, float], conf: float) -> float:
    r = conf * max(min(rates), 0.05) + 0.5 * (1 - conf)
    k, n = min(pair), sum(pair)
    return max(-(n - 2 * k) * log(r / (1 - r)), log(2.0))


def update_snps_to_use_size_cluster(
    RD: RNAData, snps_to_use: Dict[int, List[int]], size_factor: int
) -> Dict[int, List[int]]:
    return { 
        min(cluster): sorted(cluster)
        for start in snps_to_use
        for cluster in clusters_size(snps_to_use[start], RD.counts, size_factor)
    }


def update_snps_to_use_rate_cluster(
    RD: RNAData, snps_to_use: Dict[int, List[int]], rate_factor: float
) -> Dict[int, List[int]]:
    return { 
        min(cluster): sorted(cluster)
        for start in snps_to_use
        for cluster in clusters_rate(snps_to_use[start], RD.counts, rate_factor)
        if len(cluster) > 0
    }


def update_snps_to_use_rate(
    RD: RNAData, snps_to_use: Dict[int, List[int]], rate_cutoff: float
) -> Dict[int, List[int]]:
    high_confidence_starts = [ 
        start
        for start in snps_to_use
        if max(RD.rates[start]) >= rate_cutoff
    ]
    return {start: snps_to_use[start] for start in high_confidence_starts}


def update_snps_to_use_coverage(
    RD: RNAData, snps_to_use: Dict[int, List[int]], coverage_cutoff: int
) -> Dict[int, List[int]]:
    new: Dict[int, List[int]] = {}
    for start in snps_to_use:
        covered_snps = sorted([
            snp for snp in snps_to_use[start] if sum(RD.counts[snp]) > coverage_cutoff
        ])
        if len(covered_snps) > 1:
            new[covered_snps[0]] = covered_snps
    return new


def update_snps_to_use_cutoff(
    RD: RNAData, snps_to_use: Dict[int, List[int]], cutoff: float
) -> Dict[int, List[int]]:
    new: Dict[int, List[int]] = {}
    for start in snps_to_use:
        confident_snps = sorted([
            snp for snp in snps_to_use[start] if RD.LLsnp[snp] < cutoff
        ])
        if len(confident_snps) > 1:
            new[confident_snps[0]] = confident_snps
    return new


def update_snps_to_use_rate_dependent_cutoff(
    RD: RNAData, snps_to_use: Dict[int, List[int]], cutoff: int, conf: float
) -> Dict[int, List[int]]:
    new: Dict[int, List[int]] = {}
    for start in snps_to_use:
        confident_snps = sorted([
            snp for snp in snps_to_use[start] if new_score(RD.counts[snp], RD.rates[snp], conf) > cutoff
        ])
        if len(confident_snps) > 1:
            new[confident_snps[0]] = confident_snps
    return new
