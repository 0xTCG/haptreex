def count_genes(jG):
    #comp_min to genomic region
    Gdict = {}
    jGdict = {}
    
    #genomic_region to list of comp_mins mapping to it
    jGcounts = {}
    Gcounts = {}
    G = jG.G

    finals = {'g':[],'s':[],'n':[]}
    for m in G.comp_mins:
        Gdict[m] = G.data.babyRNA.SNP_to_genomic_region[m]
    for m in jG.comp_mins:
        jGdict[m] = G.data.babyRNA.SNP_to_genomic_region[m]

    for v in jGdict.values():
        jGcounts[v] = []
        Gcounts[v] = []
        
    for k in Gdict:
        Gcounts[Gdict[k]].append(k)
    for k in jGdict:
        jGcounts[jGdict[k]].append(k)

    for m in jGcounts:
        if len(Gcounts[m]) == 0:
            finals['n'].append(m)
        if len(jGcounts[m]) > len(Gcounts[m]):
            finals['g'].append(m)
        else:
            finals['s'].append(m)
    return Gdict,jGdict,Gcounts,jGcounts,finals
        
    
