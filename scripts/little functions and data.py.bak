def covered_by_at_least(D,lower):
	t=0
	for k in D:
		if k >= lower:
			t +=D[k]
	return t

def make_BC_ds(RD):
    BC = {}
    ds = {}
    for r in RD.short_read_list.values():
	k = r.keys[0]
	if BC.has_key(k):
		BC[k] +=r.count
	else:
		BC[k] = r.count
    for v in BC.values():
        if ds.has_key(v):
            ds[v]+=1
        else:
            ds[v] = 1
    return BC,ds


def coverage_of_snps_to_use(BC,STU, lower):
    t = 0
    for snps in STU.values():
        for snp in snps:
            if BC[snp]>= lower:
                t+=1
    return t
        
