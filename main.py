from common import V, DOT
from mytime import timing
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
# import chair
from getopt import gnu_getopt as getopt, GetoptError
from typing import Dict, List, Set, Tuple

print('HapTree-X v1.0 [Seq]')
random.seed(51)

def parse_args() -> Tuple[str, str, str, str, str, str]:
    try:
        if sys.argv[1] == "chair":
            opts, args = getopt(sys.argv[2:], "hv:s:t:TPC",
                ["help", "vcf=", "sam=", "threshold=", "TRANS=", "PCR=", "CONFLICT="]
            )
            vcf, sam, threshold = "", "", 1
            trans, pcr, conflict = True, False, False
            for o, a in opts:
                if o in ["-h", "--help"]: sys.exit(0)
                elif o in ['-v', '--vcf']: vcf = a
                elif o in ['-s', '--sam']: sam = a
                elif o in ['-t', '--threshold']: threshold = int(a)
                elif o in ['-T', '--TRANS']: trans = False
                elif o in ['-P', '--PCR']: pcr = True
                elif o in ['-C', '--CONFLICT']: conflict = True
                else: assert False, f"unhandled option {o}"
            if len(args) != 1:
                raise ValueError('Output argument missing')
            outputname = args[0]
            # with timing('Chair'):
            #     chair.chair(vcf, sam, outputname, threshold, trans, pcr, conflict)
            sys.exit(0)
        elif sys.argv[1] == "phase":
            opts, args = getopt(sys.argv[2:], "hv:d:r:g:i:",
                ["help", "vcf=", "dna=", "rna=", "gene=", "isoform="]
            )
            vcf, DNAfragmat, RNAfragmat, gene_data, isoforms = "", "", "", "", ""
            for o, a in opts:
                if o in ["-h", "--help"]: sys.exit(0)
                elif o in ['-v', '--vcf']: vcf = a
                elif o in ['-d', '--dna']: DNAfragmat = a
                elif o in ['-r', '--rna']: RNAfragmat = a
                elif o in ['-g', '--gene']: gene_data = a
                elif o in ['-i', '--isoform']: isoforms = a
                else: assert False, f"unhandled option {o}"
            if len(args) != 1:
                raise ValueError('Output argument missing')
            outputname = args[0]
            if vcf == "":
                raise ValueError("Error: vcf required")
            return vcf, DNAfragmat, RNAfragmat, gene_data, outputname, isoforms
        else:
            raise ValueError(f'Mode not supported')
    except ValueError as err:
        print(f'{err}')
        sys.exit(2)
    except GetoptError as err:
        print(f'{err}')
        sys.exit(2)


vcf, DNAfragmat, RNAfragmat, gene_data, outputname, isoforms = parse_args()
pair_thresh = 0.7


g = datagen.load_dna_data(vcf, [DNAfragmat], error=0.02)
GX = alg.RNA_phase(g,
    threshold=0.001,
    pair_thresh,
    error=0.02)
    #G.read_dict, G.comp_mins, G.components)

output.make_solution(GX, G, outputname)


# with timing('Main'):
#     if DNAfragmat == "" and RNAfragmat == "": # 1, 2
#         print("Need either RNA or DNA data to phase")
#     elif DNAfragmat == "" and gene_data == "": # 3
#         print("Running RNA phasing without DASE (no gene data provided)")
#         D = datagen.make_data_from_fragmat([RNAfragmat], vcf, 0.02)
#         G = graph.Graph(D)
#         GX = alg.RNA_phase(0.001, pair_thresh, 0.02, G.read_dict, G.comp_mins, G.components)
#         output.make_solution(GX, G, outputname)
#     elif DNAfragmat != "" and (RNAfragmat == "" or gene_data == ""): # 5, 6, 7
#         print("Running DNA phasing without DASE (either no RNA or no gene data provided)")
#         D: graph.Data = None
#         with timing('[1] Making fragmat'):
#             DNAfragmats = [DNAfragmat, RNAfragmat] if RNAfragmat != "" else [DNAfragmat]
#             D = datagen.make_data_from_fragmat(DNAfragmats, vcf, 0.02)
#         G: graph.Graph = None
#         with timing('[2] Making graph'):
#             G = graph.Graph(D)
#         with timing('[3] Phasing'):
#             GX = alg.RNA_phase(0.001, pair_thresh, 0.02, G.read_dict, G.comp_mins, G.components)
#             output.make_solution(GX, G, outputname)
    # else: # 4, 8
    #     RD: rna.RNAData = None
    #     with timing('[1] Making RNA fragmat'):
    #         RD = datagen.make_RNA_data_from_fragmat(gene_data, [RNAfragmat], vcf, 0.02, isoforms)
    #     with timing('[2] Stats'):
    #         stats.stats(RD, 2, 0.2, 0.6, 0, 2, 0.001, 0.2)
    #     jG: joint.JointGraph = None
    #     with timing('[3] Graph'):
    #         D = None
    #         if DNAfragmat == "":
    #             D = rna.make_data_from_RNA_data(RD)
    #         else:
    #             D = datagen.make_data_from_fragmat([DNAfragmat], vcf, 0.02, RD.long_read_list)
    #         G = graph.Graph(D)
    #         jG = joint.JointGraph(RD, G)
    #     with timing('[4] Phasing'):
    #         jGX = alg.RNA_phase(
    #             0.001, pair_thresh, 0.02, jG.read_dict, jG.comp_mins, jG.components
    #         )
    #         output.make_solution_simple(
    #             jGX,
    #             outputname,
    #             vcfChroms,
    #             vcfPositions
    #         )  # Phasing using both 1-reads and 2-reads from DNAfragmat + RNAfragmat
