import gene_class
import read_class
import rna_class
import basic_class
import global_vars
import time
import string
import re

###########################################################################
#GENES
def determine_genes(gene_data,chroms):
        print 'Loading and formatting genes'
        f = open(gene_data,'r')
        a = list(f.readlines())
        f.close()
        genes = {}
        c=0
        for x in a:
                y=x.split('\t')
                #to only include genes in chroms from VCF
                if y[0] in chroms:
                        starts = map(int,y[11].split(',')[:-1])
                        starts = [i+int(y[1])+1 for i in starts] #Convert 0-based start positions relative to gene to 1-based genomic positions
                        lengths = map(int,y[10].split(',')[:-1])
                        g = gene_class.gene(y[3],y[0],int(y[1])+1,int(y[2]),y[5],int(y[9]),lengths,starts,'','')
                        genes[c] = g
                        g.index = c
                        c+=1
        return genes

def determine_genes_gtf(gene_data,chroms):
        print 'Loading and formatting genes'
        f = open(gene_data,'r')
        a = list(f.readlines())
        #print 'done reading lines'
        f.close()
        genes = {}
        i = 0
        while a[i][0] == '#':
                i+=1
        #print 'found i'
        c = a[i:]
        #print 'c defined'
        t = time.time()
        b = {}
        for j in range(len(a[i:])):
                b[j] = tuple(c[j][:-2].split('\t'))
        #print 'done splitting'
        #print time.time() - t
        #print len(b)
        trans_list = []
        trans_dict = {}
        for j in range(len(b)):
                y = b[j]
                if y[2] == 'transcript':
                        trans_list.append(j)
                        trans_dict[j] = {'exon_starts':[], 'exon_lengths':[],'exon_count':0}

        #print len(trans_list)
        #print trans_list[:10]

        for k in range(len(trans_list)-1):
                start = trans_list[k]
                end = trans_list[k+1]
                for j in range(start,end):
                        y = b[j]
                        if y[2] == 'exon':
                                exon_start = int(y[3])
                                exon_end = int(y[4])
                                exon_length = exon_end - exon_start + 1
                                trans_dict[start]['exon_starts'].append(exon_start)
                                trans_dict[start]['exon_lengths'].append(exon_length)
                                trans_dict[start]['exon_count']+=1

        start = trans_list[-1]
        end = len(b)

        for j in range(start,end):
                y = b[j]
                if y[2] == 'exon':
                        exon_start = int(y[3])
                        exon_end = int(y[4])
                        exon_length = exon_end - exon_start + 1
                        trans_dict[start]['exon_starts'].append(exon_start)
                        trans_dict[start]['exon_lengths'].append(exon_length)
                        trans_dict[start]['exon_count']+=1

        transcripts = {}
        for j in trans_dict:
                y = b[j]
                chrom = y[0]
                if chrom in chroms:
                        bp_start = int(y[3])
                        bp_end = int(y[4])
                        sign = y[6]
                        z = y[8]
                        w = z.split('; ')
                        v = map(lambda xx:xx.split(' '),w)
                        gene_id = v[0][1][1:-1]
                        transcript_id = v[1][1][1:-1]
                        gene_type = v[2][1][1:-1]
                        exon_starts = trans_dict[j]['exon_starts']
                        exon_count = trans_dict[j]['exon_count']
                        exon_lengths = trans_dict[j]['exon_lengths']
                        transcripts[j] = gene_class.gene(transcript_id,chrom,bp_start,bp_end,sign,exon_count,
                                              exon_lengths,exon_starts,gene_id,gene_type)
                        transcripts[j].index = j

        return transcripts

def build_isodict(isoforms):
        isodict = {}
        f = open(isoforms)
        temp = map(lambda x: x.split('\t'),f.readlines())
        for i in range(1,len(temp)):
            isodict[temp[i][0]]=[float(temp[i][9]),temp[i][3]]
        return isodict


def filter_transcripts(genes,isodict):
        gene_cov = {}
        for i in isodict:
                gene_cov[isodict[i][1]]=0
        for i in isodict:
                gene_cov[isodict[i][1]]+=isodict[i][0]

        #print len(genes)
        new_genes = {}
        types = set()
        for j in genes:
                g = genes[j]
                tid = g.transcript_id
                types.add(genes[j].gene_type)
                if isodict.has_key(tid):
                        val = isodict[tid][0]
                        if val>0:
                                if val/gene_cov[isodict[tid][1]]>.01:
                                        new_genes[j] = g
        #print len(types)
        #if len(types) <100:
        #        print set(types)
        #print len(new_genes)
        return new_genes

###########################################################################
#READS

def make_read_of_frag(frag0):
        frag = frag0.split(' ')[2:]
        #takes line from fragment matrix and makes tuple formatted read
        read = ()
        qual = frag[-1][0]
        frag = frag[:-1]
        if len(frag)%2 ==1:
                print 'fragment file error'
                print frag

        for i in range(0,len(frag),2):
                key = int(frag[i])
                for char in frag[i+1]:
                        read = read + ((key-1,int(char)),)
                        key +=1
        if len(read) == 1:
                qual = ord(qual)
        else:
                qual = 1000
        return read,qual

def make_readlist_from_fragmat(fragmats):
        #translates fragment matrix into a list of reads
        print 'Loading and formatting fragments'
        a = []
        for fragmat in fragmats:
            f = open(fragmat,'r')
            A = list(f.readlines())
            f.close()
            a = a+A

        F = len(a)
        #print str(F)+' fragments'
        read_list_list = []
        i = 0
        for r in a:
                i+=1
                read,qual =make_read_of_frag(r)
                if qual >= global_vars.quality_cutoff:
                        read_list_list.append(read)
                #if not i%1000000:
                #        print i,

        print str(len(read_list_list))+ ' reads of sufficient quality'
        read_list = {}
        read_counter = {}
        i =0
        for tup_read in read_list_list:
            i+=1
            if read_counter.has_key(tup_read):
                read_counter[tup_read]+=1
            else:
                read_counter[tup_read] = 1
            #if not i%100000:
            #    print i,
        i = 0
        print str(len(read_counter)) + ' distinct reads'
        for tup in read_counter:
            read = {}
            for k,v in tup:
                read[k] = v
            read_list[i] = read_class.READ(read,read_counter[tup],i)
            #if not i%100000:
            #        print i

            i+=1
        return read_list


def make_readlist_from_fragmat_skip_1_reads(fragmats):
        #translates fragment matrix into a list of reads
        print 'Loading and formatting fragments'
        a = []
        for fragmat in fragmats:
            f = open(fragmat,'r')
            A = list(f.readlines())
            f.close()
            a = a+A

        F = len(a)
        #print str(F)+' fragments'
        read_list_list = []
        i = 0
        for r in a:
                i+=1
                read,qual =make_read_of_frag(r)
                if qual >= global_vars.quality_cutoff and len(read) > 1:
                        read_list_list.append(read)
                #if not i%1000000:
                #        print i,

        print str(len(read_list_list))+ ' reads of sufficient quality'
        read_list = {}
        read_counter = {}
        i =0
        for tup_read in read_list_list:
            i+=1
            if read_counter.has_key(tup_read):
                read_counter[tup_read]+=1
            else:
                read_counter[tup_read] = 1
            #if not i%100000:
            #    print i,
        i = 0
        print str(len(read_counter)) + ' distinct reads'
        for tup in read_counter:
            read = {}
            for k,v in tup:
                read[k] = v
            read_list[i] = read_class.READ(read,read_counter[tup],i)
            #if not i%100000:
            #        print i

            i+=1
        return read_list



###########################################################################
#make RNA_data and DNA_data objects

def positions_names_states(vcf):
        #reading data from VCF file and formatting as lists
        f = open(vcf,'r')
        a = f.readlines()
        i = 0
        while a[i][0] == '#':
                i+=1
        f.close()
        b = map(lambda x:x.split('\t'),a[i:])
        positions = {}
        states = {}
        names = {}
        chroms = {}
        for i in xrange(len(b)):
                line = b[i]
                positions[i] = int(line[1])
                names[i] = line[2] ##change when fixed files
                chroms[i] = line[0]
                state= line[9]
                states[i] = 1#sum(s)
        k = 2
        return (states,names,chroms,positions,k)


def make_RNA_data_from_fragmat(gene_data,fragmats,vcf,error,isoforms):
        ##RNA DATA
        read_list = make_readlist_from_fragmat(fragmats)
        print 'Loading VCF file'
        S,names,chroms,positions,k = positions_names_states(vcf)
        chrom_set = set(chroms.values())
        n = len(S)
        #print str(n)+ ' SNPs in VCF file'
        print 'Preparing data for ReadGraph'
        genes = determine_genes_gtf(gene_data,chrom_set)

        isodict = None
        filtered_genes = genes
        if not (isoforms == None):
                print 'Building IsoDict'
                isodict = build_isodict(isoforms)
                print 'Filtering Transcripts'
                filtered_genes = filter_transcripts(genes,isodict)

        RNA_obj = rna_class.RNA_DATA(S,genes,filtered_genes,error,read_list,positions,names,chroms,isodict)
        return RNA_obj


def make_data_from_fragmat(fragmat,vcf,error,RNA_readlist = []):
        ##regular DNA fragmat
        read_list = make_readlist_from_fragmat_skip_1_reads(fragmat)
        if len(RNA_readlist) > 0:
                max_key = max(read_list.keys())
                for zz in RNA_readlist:
                        read_list[zz+max_key] = RNA_readlist[zz]
        for r in read_list.values():
            r.special_key = r.keys[1]
            r.rates = {0:.5,1:.5}
        print 'Loading VCF file'
        S,names,chroms,positions,k = positions_names_states(vcf)
        n = len(S)
        #print str(n)+ ' SNPs in VCF file'
        print 'Preparing data for ReadGraph'
        D = basic_class.edges_from_readlist(read_list)
        data_obj = basic_class.DATA(D,S,k,error,read_list,positions,names,chroms)
        return data_obj






