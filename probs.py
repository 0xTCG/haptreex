# [Seq] DONE

from math import log
from global_vars import CONFIDENCE
from read_class import Read


def read_val_tail(
    partial_phase: dict[int, dict[int, int]],
    p: float,
    error: float,
    read_dict: dict[int, list[Read]],
    m: int,
    m_prev: int,
) -> float:
    ## we extend a phasing and then calculating how likely that extension is
    ## this looks at the prob of a particular phasing generating the set of reads which cover "end"
    ## and for which "end" is not the smallest SNP in the read. (relevant reads for the extension)
    ## this log prob gets added to the log prob of the non-extended phasing.

    end = m
    k = 2
    a = (1 - error) / (1 - (2 * error / 3.0))
    b = (error / 3.0) / (1 - (2 * error / 3.0))
    relevant_reads = read_dict[end]
    val = 0
    val2 = 0
    for read_obj in relevant_reads:
        probs = 0
        probs2 = 0
        for strand in range(k):
            pp = partial_phase[strand]
            prob = (CONFIDENCE * read_obj.rates[strand]) + (0.5 * (1 - CONFIDENCE))
            mini_read = read_obj.mini_reads[end]
            for key in mini_read:
                if mini_read[key] == pp[key]:
                    prob = prob * a
                else:
                    prob = prob * b
            if mini_read[end] == pp[end]:
                probs2 += prob / a
            else:
                probs2 += prob / b
            probs += prob
        probs = probs
        if probs == 0:
            val = -float("inf")
        else:
            val += math.log(probs) * read_obj.count
            if not end == read_obj.special_key:
                val2 += math.log(probs2) * read_obj.count

    # p is some measure of likelihood of adjacent mutations occuring together
    # since we've know extended the haplotype we need to update the prior
    if len(partial_phase[0]) > 1:
        pre_end = m_prev
        q = 1 - p
        if partial_phase[0][end] == partial_phase[0][pre_end]:
            val += math.log(p)
        else:
            val += math.log(q)

    return val - val2
