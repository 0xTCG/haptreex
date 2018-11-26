import gene_class
import basic_class
import read_class
import global_vars
import rate_finding
import time


class RNA_DATA(object):
        def __init__(self,S,genes,filtered_genes,error,read_list,positions,names,chroms,isodict):
            self.k = 2
            self.states = S
            self.initial_genes = genes
            self.genes = filtered_genes
            self.error = error
            self.all_reads = read_list
            self.positions = positions
            self.names = names
            self.chroms = chroms
            self.isodict = isodict

######################################################################

            short_read_list = {}
            long_read_list = {}
            for i in self.all_reads:
                r = self.all_reads[i]
                r.rates = {0:.5,1:.5}
                if len(r.keys)==1:
                    short_read_list[i] = r
                else:
                    long_read_list[i] = r
                    
            self.read_list = short_read_list
            self.short_read_list = short_read_list
            self.long_read_list = long_read_list

######################################################################
           
            #print set(chroms.values())
            self.nodekeys = sorted(list(self.node_keys()))
            #print len(self.nodekeys)

            d = {}
            chrom_list = {}
            for num in range(len(self.nodekeys)):
                index = self.nodekeys[num]
                d[index] = basic_class.node(self.k,index,0,[],self.positions[index],self.names[index],self.chroms[index],num)
                chrom_list[self.chroms[index]] = 0
            self.nodes = d

            self.chrom_list = sorted(chrom_list.keys())
            #print self.chrom_list
            self.phasable_positions = {chrom:[] for chrom in self.chrom_list}
            for x in self.nodekeys:
                self.phasable_positions[self.chroms[x]].append(self.positions[x])
            self.PSdict = self.make_position_dict()

######################################################################
            #print self.phasable_positions.keys()
            print 'running through genes' 
            for g in self.genes:
                self.genes[g].find_SNPs(self)

            gene_class.make_genomic_graph(self.genes)
            print 'graph made'
            
######################################################################
            #reads = self.read_list
            reads = self.all_reads
            self.comps,self.comps_dict,self.comps_mins_dict,self.SNP_to_genomic_region,self.genomic_region_to_SNPs,self.reads_by_GR =\
            gene_class.assign_reads_to_genomic_regions(self.genes,reads)
            self.reads_by_indiv_gene = self.assign_reads_to_indiv_genes()
            
######################################################################

            #Here is where we decide which types of genes to use, those with no isoforms ("no_splicing")
            #all SNPs such that if there are multiple isoforms, those SNPs fall into all of them ("final")
            
            self.final = self.dual_gene()
            self.snps_to_use = self.find_snps_to_use_all()

######################################################################

            self.reads_by_comsnps = self.find_reads_by_comsnps()

            self.read_dict = self.make_read_dict() ##1-reads only
            self.counts = self.assign_counts_to_indiv_snps()
            self.LLsnp = self.assign_LL_dif_to_snps()
            print 'assigning rates'
            self.rates  = self.assign_rates2(self.snps_to_use)[0]
            print 'rates assigned'

            current_STU0 = self.snps_to_use


######################################################################

        def node_keys(self):
            s = {}
            for r in self.read_list.values():
                for k in r.keys:
                    s[k] = 0
            return s.keys()
        
        def make_position_dict(self):
            position_to_SNP_index_dict = {chrom:{} for chrom in self.chrom_list}
            for s in self.nodekeys:
                S = self.nodes[s]
                c = S.chrom
                p = S.position
                position_to_SNP_index_dict[c][p] = S.index
            return position_to_SNP_index_dict 

######################################################################
        
 
        def assign_reads_to_indiv_genes(self):
            D = {gr:{} for gr in self.comps_mins_dict}
            for gr in self.comps_mins_dict:
                for g in self.comps[gr]:
                    D[gr][g] = []
                    snps = self.genes[g].snps
                    for r in self.reads_by_GR[gr]:
                        tf = True
                        for k in r.keys:
                            tf = (tf and (k in snps))
                        if tf:
                            D[gr][g].append(r)
            return D


        def dual_gene(self):
            D = {}
            final = {}
            for gr in self.genomic_region_to_SNPs:
                snps = self.genomic_region_to_SNPs[gr]
                start = self.comps
                D[gr] = {s:[] for s in snps}
                start = gr
                if len(self.comps[start]) == 1:
                    gene = self.genes[list(self.comps[start])[0]]
                    final[min(gene.snps)] = sorted(gene.snps)
                else:
                    for gene in self.comps[start]:
                        for snp in self.genes[gene].snps:
                            D[gr][snp].append(gene)
                            
                    for snp in D[gr]:
                        D[gr][snp] = tuple(sorted(D[gr][snp]))
                    tmp = {}
                    for s in D[gr]:
                        if not tmp.has_key(D[gr][s]):
                                tmp[D[gr][s]] = [s]
                        else:
                                tmp[D[gr][s]].append(s)
                    for snps in tmp.values():
                        final[min(snps)] = sorted(snps)
            return final
                   
                                    
                
######################################################################
### options for snps_to_use

        def find_snps_to_use_no_splicing(self):
            snps_to_use = {}
            for gr in self.genomic_region_to_SNPs:
                if len(self.comps[gr]) == 1:
                    snps = self.genomic_region_to_SNPs[gr]
                    if len(snps)>1:
                        snps_to_use[min(snps)] = snps
            return snps_to_use

        
        def find_snps_to_use_final(self):
            snps_to_use = {}
            for start in self.final:
                snps = self.final[start]
                if len(snps)>1:
                    snps_to_use[start] = snps
            return snps_to_use
                  
        def find_snps_to_use_all(self):
            snps_to_use = {}
            for start in self.genomic_region_to_SNPs:
                snps = self.genomic_region_to_SNPs[start]
                if len(snps)>1:
                    snps_to_use[min(snps)] = sorted(snps)
            return snps_to_use


######################################################################


        def find_reads_by_comsnps(self):
        # at this point we will only use reads that fall strictly within common snps
        # we should see how many long reads we arent using.
        # we should consider adding those back in to make components (earlier)
            reads_by_comsnps = {}
            for start in sorted(self.snps_to_use.keys()):
                comp = self.snps_to_use[start]
                m = start
                reads_by_comsnps[start] = []
                gr = self.SNP_to_genomic_region[start]
                for r in self.reads_by_GR[gr]:
                    tf = True
                    for k in r.keys:
                        tf = (tf and (k in comp))
                    if tf:
                        reads_by_comsnps[start].append(r)
            return reads_by_comsnps



        def assign_counts_to_indiv_snps(self):
            ###ASSUMES 1READS ONLY #1reads
            all_counts = {}
            for s in self.read_dict:
                counts = [0,0]
                for r in self.read_dict[s]:
                    counts[r.read[s]//2] += r.count
                all_counts[s] = counts
            return all_counts

        def assign_LL_dif_to_snps(self):
            LL_gr = {}
            LL = {}
            for start in self.snps_to_use:
                    for snp in self.snps_to_use[start]:
                            LL[snp] =  global_vars.score(self.counts[snp])
            return LL #LL_gr,LL

        def make_components(self,snps_to_use):
            comp_mins = []
            comps = {}
            for start in snps_to_use:
                snps = snps_to_use[start]
                if len(snps)>1:
                    for snp in snps:
                        comps[snp] = snps
                    comp_mins.append(start)

            self.components,self.comp_mins =comps,sorted(comp_mins)



        def assign_rates2(self,snps_to_use):
            times = {}
            rates = {}
            for start in snps_to_use:
                    snps = sorted(snps_to_use[start])
                    reads = self.reads_of_snps(snps)
                    t = time.time()
                    r = rate_finding.find_rates(snps,reads,.6) ##Choose adjacent snp rate
                    t = t - time.time()
                    times[t] = snps
                    for s in snps:
                        rates[s] = r
                    for read in reads:
                        read.rates = r
            self.times = times
            return rates, sorted(rates.keys())

        def make_read_dict(self):
                read_dicts = {}
                for m in self.reads_by_comsnps:
                    for r in self.reads_by_comsnps[m]:
                        if len(r.keys)==1 :
                                k = r.keys[0]
                                if read_dicts.has_key(k):
                                        read_dicts[k].append(r)
                                else:
                                        read_dicts[k] = [r]
                        else:                                        
                                for k in r.keys:
                                    new_read = read_class.READ({k:r.read[k]},r.count,-1)
                                    if read_dicts.has_key(k):
                                            read_dicts[k].append(new_read)
                                    else:
                                            read_dicts[k] = [new_read]
                return read_dicts

        def make_read_dict_full_reads(self):
                read_dict_full = {}
                for r in self.all_reads.values():
                        for k in r.keys:
                                if read_dict_full.has_key(k):
                                        read_dict_full[k].append(r)
                                else:
                                       read_dict_full[k] = [r]
                return read_dict_full

        def reads_of_snps(self,snps):
            reads = []
            for s in snps:
                #1READS
                for r in self.read_dict[s]:
                    reads.append(r)
            return reads
 


def counts_of_snps(RD,snps):
    counts  = {s:[0,0] for s in snps}
    for s in snps:
	    for r in RD.read_dict[s]:
	        counts[s][r.read[s]]+= r.count
    return counts

