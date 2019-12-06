def jGvsG_stats(jG):
    jGS = set()
    GS = set()
    for x in jG.components:
	if x in jG.G.data.babyRNA.SNP_to_genomic_region:
		jGS.add(jG.G.data.babyRNA.SNP_to_genomic_region[x])
	else:
		jGS.add(RD.SNP_to_genomic_region[x])
		
    for x in G.components:
        GS.add(G.data.babyRNA.SNP_to_genomic_region[x])
    pjG = 0
    pG = 0

    for m in jG.comp_mins:
        n = len(jG.components[m])
        pjG += (n*(n-1))/2

    for m in G.comp_mins:
        n = len(G.components[m])
        pG += (n*(n-1))/2

    jGBL = list(map(len,[jG.components[x] for x in jG.comp_mins]))
    GBL = list(map(len,[G.components[x] for x in G.comp_mins]))
    bjG = sum(map(len,[jG.components[x] for x in jG.comp_mins]))/float(len(jG.comp_mins))
    bG = sum(map(len,[G.components[x] for x in G.comp_mins]))/float(len(G.comp_mins))
    XjG = RNA_phase(jG,.001,.8)
    XG = RNA_phase(G,.001,.8)
    TjG = [switches_in_comp2(x,XjG,V,jG)[0] for x in XjG]
    TG = [switches_in_comp2(x,XG,V,G)[0] for x in XG]


    print('Total jG blocks: ', len(jG.comp_mins))
    print('Total G blocks: ', len(G.comp_mins))

    print('Breaks in jG genes: ', len(jG.comp_mins) - len(jGS))
    print('Breaks in G genes: ', len(G.comp_mins) - len(GS))

    print('SNPs phased jG: ', len(jG.components))
    print('SNPs phased G: ', len(G.components))

    print('Total edges jG: ', len(jG.components) - len(jG.comp_mins))
    print('Total edges G: ', len(G.components) - len(G.comp_mins))

    
    print('Total pairs jG: ', pjG)
    print('Total pairs G: ', pG)

    print('Mean block length jG: ', bjG)
    print('Mean block length G: ', bG)

    print('Largest blocks jG: ', sorted(jGBL)[-5:])
    print('Largest blocks G: ', sorted(GBL)[-5:])

    print('SE jG: ', sum(TjG))
    print('SE G: ', sum(TG))
