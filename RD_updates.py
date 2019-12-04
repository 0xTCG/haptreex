import math


def update_snps_to_use_size_cluster(RD,snps_to_use,size_factor):
    new = {}
    for start in snps_to_use:
        clusters0 = clusters_size(snps_to_use[start],RD.counts,size_factor)                  
        for cluster in clusters0:
            new[min(cluster)] = sorted(cluster)
    return new


def update_snps_to_use_rate_cluster(RD,snps_to_use,rate_factor):
    new = {}
    for start in snps_to_use:
        clusters0 = clusters_rate(snps_to_use[start],RD.counts,rate_factor)                  
        for cluster in clusters0:
            if len(cluster) >0:
                new[min(cluster)] = sorted(cluster)
    return new
               
def update_snps_to_use_rate(RD,snps_to_use,rate_cutoff):
    new = {}
    high_confidence_starts =[]
    for start in snps_to_use:
        h_rate = max(RD.rates[start].values())
        if h_rate>=rate_cutoff:
            high_confidence_starts.append(start)
    for start in high_confidence_starts:
            new[start] = snps_to_use[start]
    return new

def update_snps_to_use_coverage(RD,snps_to_use,coverage_cutoff):
    new = {}
    for start in snps_to_use:
        covered_snps = []
        for snp in snps_to_use[start]:
            if sum(RD.counts[snp]) > coverage_cutoff:
                covered_snps.append(snp)
        if len(covered_snps)>1:
            new[min(covered_snps)] = sorted(covered_snps)
                
    return new

def update_snps_to_use_cutoff(RD,snps_to_use,cutoff):
    new = {}
    for start in snps_to_use:
        confident_snps = []
        for snp in snps_to_use[start]:
            if RD.LLsnp[snp] < cutoff:
                confident_snps.append(snp)
        if len(confident_snps)>1:
            new[min(confident_snps)] = sorted(confident_snps)                       
    return new


def update_snps_to_use_rate_dependent_cutoff(RD,snps_to_use,cutoff,conf):
    new = {}
    for start in snps_to_use:
        confident_snps = []
        for snp in snps_to_use[start]:
            if new_score(RD.counts[snp],RD.rates[snp],conf) > cutoff:
                confident_snps.append(snp)
        if len(confident_snps)>1:
            new[min(confident_snps)] = sorted(confident_snps)                       
    return new


def update_snps_to_use_rate_dependent_cutoff_tail(RD,snps_to_use,cutoff):
    new = {}
    for start in snps_to_use:
        confident_snps = []
        for snp in snps_to_use[start]:
            if rate_tail_score(RD.counts[snp],RD.rates[snp]) > cutoff:
                confident_snps.append(snp)
        if len(confident_snps)>1:
            new[min(confident_snps)] = sorted(confident_snps)                       
    return new


def update_snps_to_use_rate_likelihood_cutoff(RD,snps_to_use,cutoff):
    new = {}
    for start in snps_to_use:
        C = [RD.counts[s] for s in snps_to_use[start]]
        m = sum(map(min,C))
        M = sum(map(max,C))
        val = score([m,M])
        if val >cutoff:
            new[start] = snps_to_use[start]                      
    return new

def update_snps_to_use_chisquare(RD,snps_to_use,chistat):
        new = {}
        for start in snps_to_use:
                C = [RD.counts[s] for s in snps_to_use[start]]
                total_counts = list(map(sum,C))
                val = chisquare(total_counts)
                if val <= chistat:
                        new[start] = snps_to_use[start]
        return new




            
def clusters_size(snps,counts,sum_factor):
    tmp_dict = {}
    for s in snps:
        tmp_dict[s] = sum(counts[s])

    order = sorted(tmp_dict,cmp = lambda x,y:cmp(tmp_dict[x],tmp_dict[y]))
    clusters = [[order[0]]]
    i = 0
    for j in range(len(order)-1):
        s0 = order[j]
        s1 = order[j+1]
        if tmp_dict[s1] < sum_factor * tmp_dict[s0]:
            clusters[i].append(s1)
        else:
            i+=1
            clusters.append([s1])
            
    return clusters

def clusters_rate(snps,counts,rate_factor):
    tmp_dict = {}
    for s in snps:
        tmp_dict[s] = max(counts[s])/float(sum(counts[s])) - .5

    order = sorted(tmp_dict,cmp = lambda x,y:cmp(tmp_dict[x],tmp_dict[y]))
    clusters = [[order[0]]]
    i = 0
    for j in range(len(order)-1):
        s0 = order[j]
        s1 = order[j+1]
        if tmp_dict[s1] < (rate_factor) + tmp_dict[s0]:
            clusters[i].append(s1)
        else:
            i+=1
            clusters.append([s1])
            
    return clusters        
            
def new_score(pair,rates,conf):
        r = min(rates.values())
        r = max(r,.05)
        r = conf*r + .5*(1-conf)
        k = min(pair)
        n  = sum(pair)
        return max((-1)*(n-2*k)*math.log(r/(1-r)),math.log(2))
