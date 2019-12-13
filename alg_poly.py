from math import log
from read import Read
from alg import prune
from common import CONFIDENCE
from time import timing
import sys

def read_val_tail(
    genes: dict[int, dict[int, int]], 
    k: int,
    p: float,
    error: float,
    relevant_reads: list[Read],
) -> float:
    tempkeys = sorted(genes[0])
    pre_end, end = tempkeys[-2], tempkeys[-1]
    a = (1 - error) / (1 - (2 * error / 3.0))
    b = (error / 3.0) / (1 - (2 * error / 3.0))
    
    val = 0.0
    for read_obj in relevant_reads: # G.read_dict[end]
        probs = 0.0
        for strand in range(k):
            prob = 1.0 # rates[strand]
            mini_read = read_obj.mini_reads[end]
            for key in mini_read:
                prob *= (a if mini_read[key] == genes[strand][key] else b)
            probs += prob
        val += log(probs / k)
    # p is supposed to be some measure of likelihood of adjacent mutations 
    # occuring together since we've know extended the haplotype 
    # we need to update the prior
    val += log(p if genes[0][end] == genes[0][pre_end] else (1 - p))
    return val


def read_val_tail_correction(
    genes: Dict[int, Dict[int, int]], 
    error: float,
    G: easy_graph
) -> float:
    # abs of log of probability of phasing (genes) given reads covering genes
    # and containing the final SNP: end, with final SNP data removed from each read
    m = len(genes[0])
    if m <= 2:
        return 0.0

    end = max(genes[0])
    a = (1 - error) / (1 - (2 * error / 3.0))
    b = (error / 3.0) / (1 - (2 * error / 3.0))
    val = 0.0
    for read_obj in G.read_dict[end]:
        if len(read_obj) <= 2:
            continue
        probs = 0.0
        for strand in range(k):
            prob = 1.0
            mini_read = read_obj.mini_reads[end]
            for key in mini_read:
                if not key == end:
                    prob *= (a if mini_read[key] == genes[strand][key] else b)
            probs += prob
        val += log(probs / k)
    return val


def branch(
    m: int, 
    LISTS: list[Dict[int, Dict[int, int]]], 
    DICT: Any, 
    G: easy_graph, 
    threshold: float, 
    p: float
) -> Any:
    # was k=2 only

    # branch all solutions to highly probable solutions on one additional SNP
    # lists include all current solutions
    # dict has solutions and their likelihhods
    if G.phase_status[m]:
        columns = [G.states[m]]
    else:
        columns = G.state_enum[G.states[m]]
    new_LISTS = []
    new_DICT = {}
    for l in LISTS:
        # extend solution l
        tup_l = tup_of_dict(l)
        costs = {c: 0 for c in columns}
        corrections = {c: 0 for c in columns}
        for c in columns:
            # adding possible orderings of alleles for new SNP
            ###this is made for diploid!!!
            temp_l = add_column_dict(m, l, c)
            costs[c] = compute_probabilities.read_val_tail(temp_l, G, p)
            corrections[c] = compute_probabilities.read_val_tail_correction(temp_l, G)
        max_cost = max(costs.values())
        new_costs = {key: costs[key] - max_cost for key in list(costs.keys())}

        for c in list(costs.keys()):
            if new_costs[c] >= math.log(threshold / float(1 - threshold)):
                # adds extension with allele ordering c if sufficiently likely
                new_l = add_column_dict(m, l, c)
                tup_new_l = tup_of_dict(new_l)
                if tup_new_l not in new_DICT:
                    # adds new solution to the dict with its new likelihood
                    ##for some dumb reason all the entries are positive instead of negative
                    new_DICT[tup_new_l] = DICT[tup_l] - costs[c] + corrections[c]
                    new_LISTS.append(new_l)
    return new_LISTS, new_DICT


def solution(start: int, G: easy_graph, threshold: float, p: float) -> List[Dict[int, Dict[int, int]]]:
    # find phasing solution for connected component starting at start
    if G.phase_status[start]:
        phase_options = [G.states[start]]
    else:
        phase_options = G.state_enum[G.states[start]]
    LISTS = [
        {i: {start: x[i]} for i in range(G.k)} for x in list(map(list, phase_options))
    ]
    DICT = dict.fromkeys(list(map(tup_of_dict, LISTS)), 0)
    for m in G.components[start][1::]:
        LISTS, DICT = branch(m, LISTS, DICT, G, threshold, p)
        LISTS, DICT = prune(m, LISTS, DICT, G)
    return LISTS

def make_genes(G: easy_graph, threshold: float, p: float) -> Tuple[Dict[int, Dict[int, Dict[int, int]]], DUMMY_NAME]:
    # phase all components of genome
    # returns those solutions and related times
    print("Phasing")
    print("Phasing completeness: ", end=" ")
    sys.stdout.flush()
    gene_dict = {}
    times = {}
    t = 0
    L = G.comp_mins[-1]
    dif = L / 10
    perc_dif = 10
    perc = 10
    temp = L / 10

    for start in G.comp_mins:
        if start > temp:
            temp += dif
            print(str(perc) + "% ", end=" ")
            sys.stdout.flush()
            perc += perc_dif
        gene_dict[start] = solution(start, G, threshold, p)[0]

    print(" ")
    return gene_dict, times

