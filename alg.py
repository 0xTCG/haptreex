import probs
import sys
import copy
import math
import random

def branch(m,m_prev,LISTS,DICT,threshold,p,error,read_dict):
    #k=2 only
        #branch all solutions to highly probable solutions on one additional SNP
        #lists include all current solutions
        #dict has list index of solutions and their likelihhods
        columns = [(0, 1), (1, 0)]
        new_LISTS = []
        new_DICT = {}
        count_list_index = 0 #counter for length of new_LISTS
        for i in range(len(LISTS)):
                partial_phase = LISTS[i]
                #extend solution
                costs = {c:0 for c in columns}
                corrections = {c:0 for c in columns}

                for c in columns:
                        #adding possible orderings of alleles for new SNP
                        #for diploid
                        extended_partial_phase= temp_add_column_dict(m,partial_phase,c)
                        costs[c] = probs.read_val_tail(partial_phase, p, error, read_dict, m, m_prev)

                max_cost = max(costs.values())
                new_costs = {key: costs[key]-max_cost for key in costs.keys()}
                used = False
                for c in costs.keys():
                        if new_costs[c] >= math.log(threshold/float(1-threshold)):
                                #adds extension with allele ordering c if sufficiently likely
                                if used: #can use the old dict once
                                    likely_extended_partial_phase = add_column_dict(m,partial_phase,c)
                                else: #need to make a copy of the dict
                                    likely_extended_partial_phase = temp_add_column_dict(m,partial_phase,c)
                                    used = True
                                #adds new solution to the dict with its new likelihood
                                #all the entries are positive instead of negative
                                new_DICT[count_list_index] = DICT[i] - costs[c]
                                count_list_index += 1
                                new_LISTS.append(likely_extended_partial_phase)
        return new_LISTS,new_DICT

def prune(m,LISTS,DICT,components):
        #prune solutions based on number of solutions and thresholds
        global counter2
        n = max(components[m])
        L = len(LISTS)
        count_list_index = 0
        new_lists = []
        new_dict = {}
        go = False

        if m == n:
                threshold = 1
                go = True
        elif L>1000:
                threshold = .1
        elif L>500:
                threshold = .05
        elif L>100:
                threshold = .01
        else:
                threshold = .001

        threshold = abs(math.log(threshold))
        min_val = min(DICT.values())
        for i in range(len(LISTS)):
                l = LISTS[i]
                if abs(DICT[i] -min_val) < threshold +.0001:
                        new_lists.append(l)
                        new_dict[count_list_index] = DICT[i]
                        count_list_index += 1
        if go:
                #output any of the most likely solutions
                    #to print all of the most likely solutions:
                #for thing in new_lists:
                #    print(tup_of_dict(thing))
                #print('$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
                new_index = random.randint(0, len(new_lists)-1)
                new_lists = [new_lists[new_index]]
                new_dict = {0:new_dict[new_index]}
        return new_lists,new_dict


def RNA_solution(start,threshold,p,error,read_dict,components):
        #find phasing solution for connected component starting at start
        init = {0:{},1:{}}
        LISTS = [init]
        DICT = {0:0}
        skipped = False #set to True to allow allele permutations
        m_prev = None
        for m in components[start]:
            if skipped:
                LISTS,DICT = branch(m,m_prev,LISTS,DICT,threshold,p,error,read_dict)
                LISTS,DICT = prune(m,LISTS,DICT,components)
            else: # specify beginning to remove allele permutations
                init = {0:{m:0},1:{m:1}}
                LISTS = [init]
                DICT = {0:0}
                skipped = True
            m_prev = m
        return LISTS

def RNA_phase(threshold,p,error,read_dict,comp_mins,components):
    phases  = {}
    for start in comp_mins:
        phases[start] = RNA_solution(start,threshold,p,error,read_dict,components)[0]
    return phases


###################################################
###### helper functions

def mycopy(phase):
    #provides a deep copy of a gene/phase
    return {0:{key:phase[0][key] for key in phase[0]}, 1:{key:phase[1][key] for key in phase[1]}}

def tup_of_dict(genes):
        #formats dictionary as a tuple so it can be used as a key
    s1 = genes[0]
    s2 = genes[1]
    keys = sorted(s1.keys())
    S1 = []
    S2 = []
    for key in keys:
        S1.append(s1[key])
        S2.append(s2[key])
    return (tuple(S1),tuple(S2))


def add_column_dict(m,l,column):
        #adds column to solution l
    new = mycopy(l)
    for i in range(len(column)):
        new[i][m] = column[i]
    return new

def temp_add_column_dict(m,l,column):
        #adds column to solution l (on exisiting solution)
    new = l
    for i in range(len(column)):
        new[i][m] = column[i]
    return new
