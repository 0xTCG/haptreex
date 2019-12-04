from pyannotate_runtime import collect_types

import global_vars
import datagen
import joint_class
import rna_class
import basic_class
import alg
import stats
import output
import prescoring
import os
import time
import sys
import argparse
import subprocess

# HapTree-X V1.0

parser = argparse.ArgumentParser(description='HapTree-X!')
parser.add_argument('--RNAfragmat',metavar='filename',type=str,help='RNAfragmat',default=None)
parser.add_argument('--DNAfragmat',metavar='filename',type=str,help='DNAfragmat',default=None)
parser.add_argument('--isoforms',metavar='filename',type=str,help='Isoform file',default=None)
parser.add_argument('--vcf',metavar='filename',type=str,help='VCF file')
parser.add_argument('--outputfolder',metavar='foldername',type=str,help='Output folder name',default='HapTreeXResults')
parser.add_argument('--gene_info',metavar='filename',type=str,help='Gene model file with GTF format',default=None)

args = parser.parse_args()

RNAfragmat = args.RNAfragmat
DNAfragmat = args.DNAfragmat
isoforms = args.isoforms
gene_data = args.gene_info
outputname = args.outputfolder
vcf = args.vcf

time_start = time.time()

if vcf == None:
    sys.exit('Error: vcf required')

if RNAfragmat == None and DNAfragmat == None:
    sys.exit('Error: at least one of RNAfragmat or DNAfragmat should be given. You can use Chair tool to generate them from SAM files.')

RNAfragmats = [RNAfragmat]
DNAfragmats = [DNAfragmat]

if not os.path.exists(outputname):
    subprocess.call(['mkdir',outputname])
pair_thresh = .7
#print 'phasing' + DNAfragmats[0]
print('Loading VCF file for scoring')
global_vars.V, global_vars.vcfChroms, global_vars.vcfPositions = prescoring.make_golden_from_true2(vcf)

if RNAfragmat == None or gene_data == None:
    print('Missing RNAfragmat, running phasing without DASE')
    if DNAfragmat == None: #Only RNA data provided without gene info
        D = datagen.make_data_from_fragmat(RNAfragmats,vcf,.02)
        G = basic_class.easy_graph(D)
        GX = alg.RNA_phase(.001,pair_thresh,.02,G.read_dict, G.comp_mins, G.components)
        output.make_solution(GX,G,outputname,'HapTreeX_noDASE_output.txt') #Running regular Haptree on 2-reads in the genome
    else:
        if RNAfragmat != None:
            DNAfragmats = [DNAfragmat, RNAfragmat]
        D = datagen.make_data_from_fragmat(DNAfragmats,vcf,.02)
        time1 = time.time()
        print(('Pt 1 took ', time1-time_start))
        G = basic_class.easy_graph(D)
        time2 = time.time()
        print(('Pt 2 took ', time2-time1))
        GX = alg.RNA_phase(.001,pair_thresh,.02,G.read_dict, G.comp_mins, G.components)
        time3 = time.time()
        print(('Pt 3 took ', time3-time2))
        output.make_solution(GX,G,outputname,'HapTreeX_noDASE_output.txt') #Running regular Haptree on 2-reads in the genome
        time4 = time.time()
        print(('Pt 4 took ', time4-time3))
        print(('Total took ', time4-time_start))
else:
    RD = datagen.make_RNA_data_from_fragmat(gene_data,RNAfragmats,vcf,.02,isoforms)
    time1 = time.time()
    print(('Pt 1 took ', time1-time_start))
    stats.stats(RD,2,.2,.6,0,2,.001,.2)
    times = time.time()
    print(('stats took ', times-time1))
    if DNAfragmat == None:
        D = basic_class.make_data_from_RNA_data(RD)
    else:
        D = datagen.make_data_from_fragmat(DNAfragmats,vcf,.02,RD.long_read_list)
    G = basic_class.easy_graph(D)
    #RDX = alg.RNA_phase(.001,pair_thresh,.02,RD.read_dict, RD.comp_mins, RD.components) #pair_thresh is prob of 00 or 11 as opposed to 01 or 10
    #output.make_solution_RNA(RDX,RD,outputname,'HapTreeX_onlyDASE_output.txt') #Phasing using only 1-reads in RNAfragmat
    #GX = alg.RNA_phase(.001,pair_thresh,.02,G.read_dict, G.comp_mins, G.components)
    #output.make_solution(GX,G,outputname,'HapTreeX_noDASE_output.txt') #Phasing using only 2-reads from RNAfragmat + DNAfragmat
    jG = joint_class.joint_graph(RD,G)
    time2 = time.time()
    print(('Pt 2 took ', time2-times))
    jGX = alg.RNA_phase(.001,pair_thresh,.02,jG.read_dict, jG.comp_mins, jG.components)
    time3 = time.time()
    print(('Pt 3 took ', time3-time2))
    output.make_solution_simple(jGX, outputname,'HapTreeX_withDASE_output.txt', global_vars.vcfChroms, global_vars.vcfPositions) #Phasing using both 1-reads and 2-reads from DNAfragmat + RNAfragmat
    time4 = time.time()
    print(('Pt 4 took ', time4-time3))
    print(('Total took ', time4-time_start))
