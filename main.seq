from common import V, DOT
from time import timing
import datagen
import joint
import rna
import graph
import alg
import stats
import output
import os
import sys
import random
from getopt import getopt, GetoptError
import chair

print 'HapTree-X v1.0 [Seq]'
random.seed(51)

def parse_args() -> tuple[str, str, str, str, str, str]:
    try:
        match sys.argv: # cannot match array pattern with non-array value
            case [_, "chair", ...]:
                opts, args = getopt(sys.argv[2:], "hv:s:t:TPC", 
                    ["help", "vcf=", "sam=", "threshold=", "TRANS=", "PCR=", "CONFLICT="]
                )
                vcf, sam, threshold = "", "", 1
                trans, pcr, conflict = True, False, False
                for o, a in opts:
                    if o in ("-h", "--help"):
                        sys.exit(0)
                    match o:
                        case '-v' or '--vcf': vcf = a
                        case '-s' or '--sam': sam = a
                        case '-t' or '--threshold': threshold = int(a)
                        case '-T' or '--TRANS': trans = False
                        case '-P' or '--PCR': pcr = True
                        case '-C' or '--CONFLICT': conflict = True
                        case _: assert False, f"unhandled option {o}"
                if len(args) != 1:
                    raise ValueError('Output argument missing')
                outputname = args[0]
                with timing('Chair'):
                    chair.chair(vcf, sam, outputname, threshold, trans, pcr, conflict)
                sys.exit(0)
            case [_, "phase", ...]:
                opts, args = getopt(sys.argv[2:], "hv:d:r:g:i:", 
                    ["help", "vcf=", "dna=", "rna=", "gene=", "isoform="]
                )
                vcf, DNAfragmat, RNAfragmat, gene_data, isoforms = "", "", "", "", ""
                for o, a in opts:
                    if o in ("-h", "--help"):
                        sys.exit(0)
                    match o:
                        case '-v' or '--vcf': vcf = a
                        case '-d' or '--dna': DNAfragmat = a
                        case '-r' or '--rna': RNAfragmat = a
                        case '-g' or '--gene': gene_data = a
                        case '-i' or '--isoform': isoforms = a
                        case _: assert False, f"unhandled option {o}"
                if len(args) != 1:
                    raise ValueError('Output argument missing')
                outputname = args[0]
                if vcf == "":
                    raise ValueError("Error: vcf required")
                return vcf, DNAfragmat, RNAfragmat, gene_data, outputname, isoforms
            case _: 
                raise ValueError(f'Mode not supported')
    except ValueError as err:
        print f'{err.message}'
        sys.exit(2)
    except GetoptError as err:
        print f'{err.message}'
        sys.exit(2)


vcf, DNAfragmat, RNAfragmat, gene_data, outputname, isoforms = parse_args()
pair_thresh = 0.7

def make_golden_from_true2(
    filename: str
) -> tuple[dict[int, dict[int, int]], list[str], list[str]]:
    f = open(filename, "r")
    a = f.readlines()
    f.close()
    i = 0
    while a[i][0] == "#":
        i += 1
    c = [x.split()[9][:3] for x in a[i:]]
    vChroms = [x.split()[0] for x in a[i:]]
    vPositions = [x.split()[1] for x in a[i:]]
    d = {0: dict[int, int](), 1: dict[int, int]()}
    for j in range(len(c)):
        if c[j][1] == "|":
            d[0][j] = int(c[j][0])
            d[1][j] = int(c[j][2])
        else:
            d[0][j] = DOT
            d[1][j] = DOT
    return d, vChroms, vPositions
vcfChroms, vcfPositions = list[str](), list[str]()
with timing(f"[0] Loading gold VCF {vcf}"):
    V, vcfChroms, vcfPositions = make_golden_from_true2(vcf)

with timing('Main'):
    if DNAfragmat == "" and RNAfragmat == "": # 1, 2
        print "Need either RNA or DNA data to phase"
    elif DNAfragmat == "" and gene_data == "": # 3
        print "Running RNA phasing without DASE (no gene data provided)"
        D = datagen.make_data_from_fragmat([RNAfragmat], vcf, 0.02)
        G = graph.Graph(D)
        GX = alg.RNA_phase(0.001, pair_thresh, 0.02, G.read_dict, G.comp_mins, G.components)
        output.make_solution(GX, G, outputname)
    elif DNAfragmat != "" and (RNAfragmat == "" or gene_data == ""): # 5, 6, 7
        print "Running DNA phasing without DASE (either no RNA or no gene data provided)"
        D: graph.Data = None
        with timing('[1] Making fragmat'):
            DNAfragmats = [DNAfragmat, RNAfragmat] if RNAfragmat != "" else [DNAfragmat]
            D = datagen.make_data_from_fragmat(DNAfragmats, vcf, 0.02)
        G: graph.Graph = None
        with timing('[2] Making graph'):
            G = graph.Graph(D)
        with timing('[3] Phasing'):
            GX = alg.RNA_phase(0.001, pair_thresh, 0.02, G.read_dict, G.comp_mins, G.components)
            output.make_solution(GX, G, outputname)
    else: # 4, 8
        RD: rna.RNAData = None
        with timing('[1] Making RNA fragmat'):
            RD = datagen.make_RNA_data_from_fragmat(gene_data, [RNAfragmat], vcf, 0.02, isoforms)
        with timing('[2] Stats'):
            stats.stats(RD, 2, 0.2, 0.6, 0, 2, 0.001, 0.2)
        jG: joint.JointGraph = None
        with timing('[3] Graph'):
            D = None
            if DNAfragmat == "":
                D = rna.make_data_from_RNA_data(RD)
            else:
                D = datagen.make_data_from_fragmat([DNAfragmat], vcf, 0.02, RD.long_read_list)
            G = graph.Graph(D)
            jG = joint.JointGraph(RD, G)
        with timing('[4] Phasing'):
            jGX = alg.RNA_phase(
                0.001, pair_thresh, 0.02, jG.read_dict, jG.comp_mins, jG.components
            )
            output.make_solution_simple(
                jGX,
                outputname,
                vcfChroms,
                vcfPositions
            )  # Phasing using both 1-reads and 2-reads from DNAfragmat + RNAfragmat
