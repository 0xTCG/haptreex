from math import log
from probs import read_val_tail
from read_class import Read


def add_column_dict(
    m: int, l: dict[int, dict[int, int]], column: tuple[int, int]
) -> dict[int, dict[int, int]]:
    # adds column to solution l
    new = { 0: {key: l[0][key] for key in l[0]},
            1: {key: l[1][key] for key in l[1]} }
    for i in range(len(column)):
        new[i][m] = column[i]
    return new


def temp_add_column_dict(
    m: int, l: dict[int, dict[int, int]], column: tuple[int, int]
) -> dict[int, dict[int, int]]:
    # adds column to solution l (on exisiting solution)
    new = l
    for i in range(len(column)):
        new[i][m] = column[i]
    return new


# def branch(m: int,m_prev: Optional[int],LISTS: Union[list[dict[int, dict[int, int]]], list[dict[int, DUMMY_NAME]]],DICT: dict[int, Union[int, float]],threshold: float,p: float,error: float,read_dict: dict[int, list[Read]]) -> tuple[list[dict[int, dict[int, int]]], dict[int, float]]:
def branch(
    m: int,
    m_prev: int,
    lists: list[dict[int, dict[int, int]]],
    table: dict[int, Union[int, float]],
    threshold: float,
    p: float,
    error: float,
    read_dict: dict[int, list[Read]],
) -> tuple[list[dict[int, dict[int, int]]], dict[int, float]]:
    # k=2 only
    # branch all solutions to highly probable solutions on one additional SNP
    # lists include all current solutions
    # dict has list index of solutions and their likelihhods
    columns = [(0, 1), (1, 0)]
    new_lists = list[dict[int, dict[int, int]]]()
    new_table = dict[int, float]()
    count_list_index = 0  # counter for length of new_lists
    for i in range(len(lists)):
        partial_phase = lists[i]
        # extend solution
        costs = {c: 0 for c in columns}
        corrections = {c: 0 for c in columns}

        for c in columns:
            # adding possible orderings of alleles for new SNP for diploid
            extended_partial_phase = temp_add_column_dict(m, partial_phase, c)
            costs[c] = read_val_tail(partial_phase, p, error, read_dict, m, m_prev)

        max_cost = max(costs.values())
        new_costs = {key: costs[key] - max_cost for key in list(costs.keys())}
        used = False
        for c in costs.keys():
            if new_costs[c] >= log(threshold / (1 - threshold)):
                # adds extension with allele ordering c if sufficiently likely
                if used:  # can use the old dict once
                    likely_extended_partial_phase = add_column_dict(m, partial_phase, c)
                else:  # need to make a copy of the dict
                    likely_extended_partial_phase = temp_add_column_dict(
                        m, partial_phase, c
                    )
                    used = True
                # adds new solution to the dict with its new likelihood
                # all the entries are positive instead of negative
                new_table[count_list_index] = table[i] - costs[c]
                count_list_index += 1
                new_lists.append(likely_extended_partial_phase)
    return new_lists, new_table


def prune(
    m: int,
    lists: list[dict[int, dict[int, int]]],
    table: dict[int, float],
    components: dict[int, list[int]],
) -> tuple[list[dict[int, dict[int, int]]], dict[int, float]]:
    # prune solutions based on number of solutions and thresholds
    n = max(components[m])
    L = len(lists)
    count_list_index = 0
    new_lists = list[dict[int, dict[int, int]]]()
    new_dict = dict[int, float]()
    go = False
    threshold = 1.0
    if m == n:
        threshold = 1
        go = True
    elif L > 1000:
        threshold = 0.1
    elif L > 500:
        threshold = 0.05
    elif L > 100:
        threshold = 0.01
    else:
        threshold = 0.001

    threshold = abs(log(threshold))
    min_val = min(table.values())
    for i in range(len(lists)):
        l = lists[i]
        if abs(table[i] - min_val) < threshold + 0.0001:
            new_lists.append(l)
            new_dict[count_list_index] = table[i]
            count_list_index += 1
    if go:
        # output any of the most likely solutions
        # to print all of the most likely solutions:
        # for thing in new_lists:
        #    print(tup_of_dict(thing))
        # print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
        new_index = 0  # random.randint(0, len(new_lists)-1)
        new_lists = [new_lists[new_index]]
        new_dict = {0: new_dict[new_index]}
    return new_lists, new_dict


def RNA_solution(
    start: int,
    threshold: float,
    p: float,
    error: float,
    read_dict: dict[int, list[Read]],
    components: dict[int, list[int]],
) -> list[dict[int, dict[int, int]]]:
    # find phasing solution for connected component starting at start
    init = {0: dict[int, int](), 1: dict[int, int]()}
    lists = [init]
    table = {0: 0}
    skipped = True  # set to True to allow allele permutations
    m_prev = -1 # SEQ/TODO: check! originally  None
    for m in components[start]:
        if skipped:
            lists, table = branch(m, m_prev, lists, table, threshold, p, error, read_dict)
            lists, table = prune(m, lists, table, components)
        else:  # specify beginning to remove allele permutations
            init = {0: {m: 0}, 1: {m: 1}}
            lists = [init]
            table = {0: 0}
            skipped = True
        m_prev = m
    return lists


def RNA_phase(
    threshold: float,
    p: float,
    error: float,
    read_dict: dict[int, list[Read]],
    comp_mins: list[int],
    components: dict[int, list[int]],
):
    phases = dict[int, dict[int, dict[int, int]]]()
    for s in comp_mins:
        phases[s] = RNA_solution(start, threshold, p, error, read_dict, components)[0]
    return phases
