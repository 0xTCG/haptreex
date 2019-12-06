from math import log
from rna import RNAData


def clusters_size(
    snps: list[int], counts: dict[int, list[int]], sum_factor: int
) -> list[list[int]]:
    tmp_dict = dict[int, float]()
    for s in snps:
        tmp_dict[s] = sum(counts[s])

    order = sorted(tmp_dict, key=lambda x: tmp_dict[x])
    clusters = [[order[0]]]
    i = 0
    for j in range(len(order) - 1):
        s0 = order[j]
        s1 = order[j + 1]
        if tmp_dict[s1] < sum_factor * tmp_dict[s0]:
            clusters[i].append(s1)
        else:
            i += 1
            clusters.append([s1])

    return clusters


def clusters_rate(
    snps: list[int], counts: dict[int, list[int]], rate_factor: float
) -> list[list[int]]:
    tmp_dict = dict[int, float]()
    for s in snps:
        tmp_dict[s] = max(counts[s]) / float(sum(counts[s])) - 0.5

    order = sorted(tmp_dict, key=lambda x: tmp_dict[x])
    clusters = [[order[0]]]
    i = 0
    for j in range(len(order) - 1):
        s0 = order[j]
        s1 = order[j + 1]
        if tmp_dict[s1] < (rate_factor) + tmp_dict[s0]:
            clusters[i].append(s1)
        else:
            i += 1
            clusters.append([s1])

    return clusters


def new_score(pair: list[int], rates: dict[int, float], conf: float) -> float:
    r = min(rates.values())
    r = max(r, 0.05)
    r = conf * r + 0.5 * (1 - conf)
    k = min(pair)
    n = sum(pair)
    return max((-1) * (n - 2 * k) * log(r / (1 - r)), log(2))


def update_snps_to_use_size_cluster(
    RD: RNAData, snps_to_use: dict[int, list[int]], size_factor: int
) -> dict[int, list[int]]:
    new = dict[int, list[int]]()
    for start in snps_to_use:
        clusters0 = clusters_size(snps_to_use[start], RD.counts, size_factor)
        for cluster in clusters0:
            new[min(cluster)] = sorted(cluster)
    return new


def update_snps_to_use_rate_cluster(
    RD: RNAData, snps_to_use: dict[int, list[int]], rate_factor: float
) -> dict[int, list[int]]:
    new = dict[int, list[int]]()
    for start in snps_to_use:
        clusters0 = clusters_rate(snps_to_use[start], RD.counts, rate_factor)
        for cluster in clusters0:
            if len(cluster) > 0:
                new[min(cluster)] = sorted(cluster)
    return new


def update_snps_to_use_rate(
    RD: RNAData, snps_to_use: dict[int, list[int]], rate_cutoff: float
) -> dict[int, list[int]]:
    new = dict[int, list[int]]()
    high_confidence_starts = list[int]()
    for start in snps_to_use:
        h_rate = max(RD.rates[start].values())
        if h_rate >= rate_cutoff:
            high_confidence_starts.append(start)
    for start in high_confidence_starts:
        new[start] = snps_to_use[start]
    return new


def update_snps_to_use_coverage(
    RD: RNAData, snps_to_use: dict[int, list[int]], coverage_cutoff: int
) -> dict[int, list[int]]:
    new = dict[int, list[int]]()
    for start in snps_to_use:
        confident_snps = list[int]()
        for snp in snps_to_use[start]:
            if sum(RD.counts[snp]) > coverage_cutoff:
                covered_snps.append(snp)
        if len(covered_snps) > 1:
            new[min(covered_snps)] = sorted(covered_snps)
    return new


def update_snps_to_use_cutoff(
    RD: RNAData, snps_to_use: dict[int, list[int]], cutoff: float
) -> dict[int, list[int]]:
    new = dict[int, list[int]]()
    for start in snps_to_use:
        confident_snps = list[int]()
        for snp in snps_to_use[start]:
            if RD.LLsnp[snp] < cutoff:
                confident_snps.append(snp)
        if len(confident_snps) > 1:
            new[min(confident_snps)] = sorted(confident_snps)
    return new


def update_snps_to_use_rate_dependent_cutoff(
    RD: RNAData, snps_to_use: dict[int, list[int]], cutoff: int, conf: float
) -> dict[int, list[int]]:
    new = dict[int, list[int]]()
    for start in snps_to_use:
        confident_snps = list[int]()
        for snp in snps_to_use[start]:
            if new_score(RD.counts[snp], RD.rates[snp], conf) > cutoff:
                confident_snps.append(snp)
        if len(confident_snps) > 1:
            new[min(confident_snps)] = sorted(confident_snps)
    return new
