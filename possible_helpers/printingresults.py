def what_to_print(RD,X,V):

    results = []

    #phasing info
    snps_by_chrom = {c:[] for c in RD.chrom_list}
    unknowns_by_chrom = {c:[] for c in RD.chrom_list}
    switches_by_chrom = {c:[] for c in RD.chrom_list}
    for x in X:
        c = RD.nodes[x].chrom
        snps_by_chrom[c].append(len(X[x]))
        unknowns = list(V[x][0].values()).count('.')
        unknows_by_chrom[c].append(unknowns)
        switches = switches_by_comp(x,X,V,RD)[0]
        switches_by_chrom[c].append(switches)
    
    results.append(sum(map(sum,list(snps_by_chrom.values()))))
    results.append(sum(map(sum,list(unkowns_by_chrom.values()))))
    results.append(sum(map(sum,list(switches_by_chrom.values()))))
    results.append(snps_by_chrom)
    results.append(unknowns_by_chrom)
    results.append(switches_by_chrom)


def gene_data_init(RD):
    results = []
    all_genes = len(RD.genes)
    genes_covered = len(RD.comps)
    results.append(['all genes: ', all_genes])
    results.append(['genes covered: ', genes_covered])
    genomic_regions = list(map(len,list(RD.genomic_region_to_SNPs().values())))
    results.append(['genomic regions: ',len(genomic_regions)])
    results.append(['genomic region sizes: ',sorted(genomic_regions)])
    results.append(['counted genomic region sizes: ',list(map(genomic_regions.count,list(range(1,max(genomic_regions)+1))))])
    
    results.append(['total snps seen: ',sum(genomic_regions)])
    snp_cov = []
    for snp in RD.read_dict:
        counter = 0
        if len(RD.read_dict[snp])>0:
            for r in RD.read_dict[snp]:
                counter += r.count
        snp_cov.append(counter)
    results.append(['coverage of snps: ',sorted(snp_cov)])
    results.append(['counted coverage: ',list(map(snp_cov.count,list(range(1,100))))])
    results.append(['max coverage: ',max(snp_cov)])
    

        
        

    
    

    
        
        
    
