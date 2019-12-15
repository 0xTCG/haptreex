from math import log
from read import Read
from common import CONFIDENCE
from mytime import timing
from random import randint
import sys
from copy import copy
from typing import Tuple, Dict, List, Set


def read_val_tail(
    phase: List[Dict[SNP, int]],
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

    assert len(phase) == 2 # Diploid for now

    a = (1 - error) / (1 - (2 * error / 3.0))
    b = (error / 3.0) / (1 - (2 * error / 3.0))
    val, val2 = 0.0, 0.0
    for read in relevant_reads:
        probs, probs2 = 0.0, 0.0
        for i, pp in enumerate(phase):
            prob = (CONFIDENCE * read.rates[i]) + (0.5 * (1 - CONFIDENCE))
            for snp in read.snps:
                if snp <= last_snp:
                    prob *= (a if read.snps[snp] == pp[snp] else b)
            probs2 += prob / (a if read.snps[last_snp] == pp[last_snp] else b)
            probs += prob
        if probs == 0:
            val = -float("inf")
        else:
            val += log(probs) * read.count
            if not last_snp == read.special_snp:
                val2 += log(probs2) * read.count

    # pair_threshold is some measure of likelihood of adjacent mutations occuring together
    # since we've know extended the haplotype we need to update the prior
    if len(phase[0]) > 1:
        assert prev_snp
        val += log(
            pair_threshold
            if phase[0][last_snp] == phase[0][prev_snp]
            else (1 - pair_threshold)
        )

    return val - val2


def branch(
    phase: Tuple[List[List[Dict[SNP, int]]], List[float]],
    snp: SNP,
    prev_snp: Optional[SNP],
    relevant_reads: list[Read],
    threshold: float,
    pair_threshold: float,
    error: float,
):
    """
    k=2 only
    branch all solutions to highly probable solutions on one additional SNP
    lists include all current solutions
    Dict has List index of solutions and their likelihoods
    """

    lists, table = phase
    new_lists, new_table = [], []
    configs, costs = [(0, 1), (1, 0)], [0.0, 0.0]
    for pi, partial in enumerate(lists):
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
                extended_phase = partial
                if used:
                    extended_phase = (copy(partial[0]), copy(partial[1]))
                else:
                    used = True
                extended_phase[0][snp], extended_phase[1][snp] = configs[i]
                new_table.append(table[pi] - costs[i])
                new_lists.append(extended_phase)
    return new_lists, new_table


def prune(
    phase: Tuple[List[List[Dict[SNP, int]]], List[float]],
    is_last: bool
):
    lists, table = phase
    # prune solutions based on number of solutions and thresholds
    threshold = 1.0
    if is_last: threshold = 1.0
    elif len(lists) > 1000: threshold = 0.1
    elif len(lists) > 500: threshold = 0.05
    elif len(lists) > 100: threshold = 0.01
    else: threshold = 0.001
    threshold = abs(log(threshold)) + 0.0001

    min_val = min(table)
    nlists: List[List[Dict[SNP, int]]] = []
    ntable: List[float] = []
    for i, l in enumerate(lists):
        if abs(table[i] - min_val) < threshold:
            nlists.append(l)
            ntable.append(table[i])
    if is_last:
        ni = randint(0, len(nlists) - 1)
        nlists, ntable = [nlists[ni]], [ntable[ni]]
    return nlists, ntable


def parallel(
    graph: Graph,
    component: Tuple[int, SNP],
    threshold: float,
    pair_threshold: float,
    error: float,
    # read_dict: Dict[int, List[Read]],
    # components: Dict[int, List[int]],
    output: List[Tuple[SNP, List[Dict[SNP, int]]]]
):
    idx, root = component
    # Ploidy Phase: [ SNP -> Allele ]
    phase: Tuple[List[List[Dict[SNP, int]]], List[float]] = [[{}, {}]], [0.0]
    skipped = True  # Set to True to allow allele permutations; changed from original False
    prev_snp: SNP = None # SEQ/TODO: check! originally None
    for snp in graph.components[root]:
        if skipped:
            phase = branch(phase, snp, prev_snp, graph.reads[snp], threshold, pair_threshold, error)
            phase = prune(phase, snp == graph.components[snp][-1])
        else:  # specify beginning to remove allele permutations
            phase = [[{}, {}]], [0.0]
            skipped = True
        prev_snp = snp
    output[idx] = (root, phase[0][0])


def RNA_phase(
    graph: Graph,
    threshold: float,
    pair_threshold: float,
    error: float
):
    threshold = log(threshold / (1 - threshold))
    output = [None for _ in graph.component_roots] #S [(0, (Dict[int, int](), Dict[int, int]())) for _ in comp_mins]
    print(f'Phasing {len(graph.component_roots)} components...')
    for i in enumerate(graph.component_roots):
        parallel(i, graph, threshold, pair_threshold, error)
    return dict(output) # (component) -> phase