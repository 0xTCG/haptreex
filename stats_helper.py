from math import log
from rna import RNAData
from typing import Tuple, Dict, List, Set




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
