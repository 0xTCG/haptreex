from math import log
from read import Read
from common import CONFIDENCE
from mytime import timing
from random import choice
import sys
from copy import copy
from typing import Tuple, Dict, List, Set, NamedTuple


class Phase(NamedTuple):
    haplotypes: List[Dict[SNP, int]] # A dictionary for each haplotype
    score: float


def read_val_tail(
    haplotypes: List[Dict[SNP, int]]
    pair_threshold: float,
    error: float,
    relevant_reads: List[Read],
    last_snp: SNP,
    prev_snp: Optional[SNP]
) -> float:
    """
    we extend a phasing and then calculating how likely that extension is
    this looks at the prob of a particular phasing generating the Set of reads
    which cover "end" ( = "last_snp") and for which "end" is not the smallest SNP in the read.
    (relevant reads for the extension)
    this log prob gets added to the log prob of the non-extended phasing.
    """

    assert len(haplotypes) == 2 # Diploid for now

    a = (1 - error) / (1 - (2 * error / 3.0))
    b = (error / 3.0) / (1 - (2 * error / 3.0))
    val, val2 = 0.0, 0.0
    for read in relevant_reads:
        probs, probs2 = 0.0, 0.0
        for hi, haplo in enumerate(haplotypes):
            prob = (CONFIDENCE * read.rates[hi]) + (0.5 * (1 - CONFIDENCE))
            for snp in read.snps:
                if snp <= last_snp:
                    prob *= (a if read.snps[snp] == haplo[snp] else b)
            probs2 += prob / (a if read.snps[last_snp] == haplo[last_snp] else b)
            probs += prob
        if probs == 0:
            val = -float("inf")
        else:
            val += log(probs) * read.count
            if not last_snp == read.special_snp:
                val2 += log(probs2) * read.count

    # pair_threshold is some measure of likelihood of adjacent mutations occuring together
    # since we've know extended the haplotype we need to update the prior
    if len(haplotypes[0]) > 1:
        assert prev_snp
        val += log(
            pair_threshold
            if haplotypes[0][last_snp] == haplotypes[0][prev_snp]
            else (1 - pair_threshold)
        )
    return val - val2


def branch(
    phases: List[Phase],
    snp: SNP,
    prev_snp: Optional[SNP],
    relevant_reads: list[Read],
    threshold: float,
    pair_threshold: float,
    error: float,
) -> List[Phase]:
    """
    k=2 only
    branch all solutions to highly probable solutions on one additional SNP
    lists include all current solutions
    Dict has List index of solutions and their likelihoods
    """

    new_phases = []
    configs, costs = [(0, 1), (1, 0)], [0.0, 0.0]
    for partial, score in enumerate(phases):
        assert len(partial) == 2
        # extend solution
        for i in range(2):
            # adding possible orderings of alleles for new SNP for diploid
            partial[0][snp], partial[1][snp] = configs[i]
            costs[i] = read_val_tail(partial, pair_threshold, error, relevant_reads, snp, prev_snp)
        max_cost = max(costs[0], costs[1])
        used = False
        for i in range(2):
            if costs[i] - max_cost >= threshold:
                # adds extension with allele ordering c if sufficiently likely
                extended = partial
                if used:
                    extended = (copy(partial[0]), copy(partial[1]))
                else:
                    used = True
                extended[0][snp], extended[1][snp] = configs[i]
                new_phases.append(Phase(extended, score - costs[i]))
    return new_phases


def prune(
    phases: List[Phase],
    is_last: bool
) -> List[Phase]:
    lists, table = phase
    # prune solutions based on number of solutions and thresholds
    threshold = 1.0
    if is_last: threshold = 1.0
    elif len(phases) > 1000: threshold = 0.1
    elif len(phases) > 500: threshold = 0.05
    elif len(phases) > 100: threshold = 0.01
    else: threshold = 0.001
    threshold = abs(log(threshold)) + 0.0001

    min_val = min(table)
    nphases = [p for p in phases if abs(p.score - min_val) < threshold]
    if is_last:
        nphases = choice(nphases)
    return nphases


def parallel(
    graph: Graph,
    root: SNP,
    threshold: float,
    pair_threshold: float,
    error: float,
) -> Phase:
    phases = [Phase([{}, {}], 0.0)]
    skipped = True  # Set to True to allow allele permutations; changed from original False
    prev_snp: SNP = None # SEQ/TODO: check! originally None
    for snp in graph.components[root]:
        if skipped:
            phases = branch(phases, snp, prev_snp, graph.reads[snp], threshold, pair_threshold, error)
            phases = prune(phases, snp == graph.components[snp][-1])
        else:  # specify beginning to remove allele permutations
            phases = [Phase([{}, {}], 0.0)]
            skipped = True
        prev_snp = snp
    return phases[0]


def phase(
    graph: Graph,
    threshold: float,
    pair_threshold: float,
    error: float
) -> Dict[SNP, Phase]:
    threshold = log(threshold / (1 - threshold))
    print(f'Phasing {len(graph.component_roots)} components...')
    result = {}
    for root in enumerate(graph.component_roots):
        result[r] = parallel(graph, root, threshold, pair_threshold, error)
    return result