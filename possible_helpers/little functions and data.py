def covered_by_at_least(D,lower):
	t=0
	for k in D:
		if k >= lower:
			t +=D[k]
	return t

def make_BC_ds(RD):
    BC = {}
    ds = {}
    for r in list(RD.short_read_list.values()):
	k = r.keys[0]
	if k in BC:
		BC[k] +=r.count
	else:
		BC[k] = r.count
    for v in list(BC.values()):
        if v in ds:
            ds[v]+=1
        else:
            ds[v] = 1
    return BC,ds


def coverage_of_snps_to_use(BC,STU, lower):
    t = 0
    for snps in list(STU.values()):
        for snp in snps:
            if BC[snp]>= lower:
                t+=1
    return t
        
