import bisect

class gene(object):
        def __init__(self,name,chrom,bp_start,bp_end,sign,exon_count,exon_lengths,exon_starts,gene_id,gene_type):
                self.name = name
                self.transcript_id = name
                self.gene_id = gene_id
                self.gene_type = gene_type
                self.chrom = chrom
                self.bp_start = bp_start #Assumes that constructor calls sends in the start positions 1-based (both for BED and GTF)
                self.bp_end = bp_end  
                self.exon_lengths = exon_lengths
                self.exon_starts = exon_starts
                self.sign = sign
                self.exon_count = exon_count
                self.ranges =  []
                #FOR GTF
                for i in range(exon_count):
                        self.ranges.append((self.exon_starts[i],self.exon_starts[i]+self.exon_lengths[i]-1))

        def find_SNPs(self,RD):
                
                s = self.find_next_SNP(self.bp_start,RD)

                if s ==-1:
                        snps = {i:[] for i in range(self.exon_count)}
                        self.snp_dict = snps
                        self.snps = []
                else:
                        S =RD.nodes[s]
                        num = S.num
                        ALL_SNPS = []
                        while S.position < self.bp_end:
                                ALL_SNPS.append(S)
                                num+=1
                                if num < len(RD.nodekeys):
                                    S = RD.nodes[RD.nodekeys[num]]
                                else:
                                    break
                        L = len(ALL_SNPS)
                        r_i = 0
                        s_i = 0
                       
                        snps = {i:[] for i in range(self.exon_count)}

                        while s_i <L:
                                S = ALL_SNPS[s_i]
                                r = self.ranges[r_i]
                                tf = in_range(S,r)
                                if (tf[0] and tf[1]):
                                        snps[r_i].append(S)
                                        s_i+=1
                                elif tf[0]:
                                        r_i +=1
                                elif tf[1]:
                                        s_i+=1
                        self.snp_dict = snps
                        self.snps = []
                        for x in self.snp_dict.values():
                                for y in x:
                                        self.snps.append(y.index)
                                        
                                
                                        
        def __eq__(self,other):
                return self.name == other.name
                        
        def __hash__(self):
                return hash(self.name)
                
        def find_next_SNP(self,b,RD):
                if not RD.phasable_positions.has_key(self.chrom):
                        return -1
                else:
                        positions = RD.phasable_positions[self.chrom]
                        L = len(positions)
                        PSdict = RD.PSdict[self.chrom]
                        nodekeys = RD.nodekeys
                        h=bisect.bisect_left(positions,b)

                        if h < L:
                                i = positions[h]
                                j = PSdict[i]
                                return j
                        else:
                                return -1


        

                
def in_range(S,r):
        p = S.position
        return (p >=r[0],p<r[1]) 
                




  
def gene_by_snp_dict(GG):
    s_dict = {}
    for g in GG.values():
        for s in g.snps:
            if s not in s_dict:
                s_dict[s] = {g.index}
            else:
                s_dict[s].add(g.index)
    return s_dict

def make_genomic_graph(GG): #making neighbours
    s_dict = gene_by_snp_dict(GG) #snp to gene
    for x in GG:
        g = GG[x]
        g.neighbors = set()
    for x in GG:
        g = GG[x]
        for s in g.snps:
            for y in s_dict[s]:
                g.neighbors.add(y)
 
    
    
        
def build_genomic_comp(start, GG):
        #builds connected component containing the node start
        #returns ordered list of nodes in component
    s = {start.index}
    old_batch = {start.index}
    while len(old_batch) > 0:
            new_batch = set()
            for g in old_batch:
                    new_batch = new_batch.union(GG[g].neighbors)
            new_batch = new_batch.difference(s)
            s=s.union(new_batch)
            old_batch = new_batch
    return s


def build_all_genomic_comps(GG):
        #print("GG", GG)

        ##For regular comps not genomic:
        #builds connected component dictionaries
        #comps maps each node to a list of nodes it is connected to
        #comps_dict maps each node to the smallest node it is connected to 
        comps = {}
        comps_dict = {}
        for start in GG:
            if not GG[start].snps == []:
                if not comps_dict.has_key(start):
                    s = build_genomic_comp(GG[start], GG)
                    comps[start] = s
                    m = min(s)
                    for g in s:
                        comps_dict[g] = m
                        comps[g] = s
        return comps,comps_dict



def assign_reads_to_genomic_regions(genes,reads):
    comps,comps_dict = build_all_genomic_comps(genes)
    comps_mins_dict = gene_comps_by_mins(comps,comps_dict)
    S_G,G_S = snp_to_genomic_region(comps_mins_dict, genes)
    reads_by_gene = {m:[] for m in comps_mins_dict}
    reads_by_gene['Over Seen'] = []
    reads_by_gene['Not Seen'] = []
    
    for R in reads.values():
        regions = set()
        for s in R.keys:
            if S_G.has_key(s):
                regions.add(S_G[s])
            else:
                regions.add('Not Seen')
        if len(regions) == 1:
            m = list(regions)[0]
            R.region = m
            reads_by_gene[m].append(R)
        else:
            reads_by_gene['Over Seen'].append(R)

    return comps,comps_dict,comps_mins_dict,S_G,G_S,reads_by_gene

        
def gene_comps_by_mins(comps,comps_dict):
    mins = list(set(comps_dict.values()))
    gcbm = {}
    for m in mins:
        gcbm[m]  = comps[m]
    return gcbm

def snp_to_genomic_region(gcmb, GG):
    SNP_GR = {}
    for m in gcmb:
        for g in gcmb[m]:
            for s in GG[g].snps:
                SNP_GR[s] = m
    GR_SNP = {m:[] for m in gcmb}
    for s in SNP_GR:
        GR_SNP[SNP_GR[s]].append(s)
    for m in GR_SNP:
        GR_SNP[m] = sorted(GR_SNP[m])
    return SNP_GR,GR_SNP
