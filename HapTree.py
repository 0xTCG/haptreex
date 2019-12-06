from global_vars import V
import datagen
import joint_class
import rna_class
import basic_class
import alg
import stats
import output
import os
import time
import sys
# import argparse
# import subprocess

# HapTree-X V1.0

# parser = argparse.ArgumentParser(description="HapTree-X!")
# parser.add_argument(
#     "--RNAfragmat", metavar="filename", type=str, help="RNAfragmat", default=None
# )
# parser.add_argument(
#     "--DNAfragmat", metavar="filename", type=str, help="DNAfragmat", default=None
# )
# parser.add_argument(
#     "--isoforms", metavar="filename", type=str, help="Isoform file", default=None
# )
# parser.add_argument("--vcf", metavar="filename", type=str, help="VCF file")
# parser.add_argument(
#     "--outputfolder",
#     metavar="foldername",
#     type=str,
#     help="Output folder name",
#     default="HapTreeXResults",
# )
# parser.add_argument(
#     "--gene_info",
#     metavar="filename",
#     type=str,
#     help="Gene model file with GTF format",
#     default=None,
# )

# args = parser.parse_args()

RNAfragmat = sys.argv[1] # args.RNAfragmat
DNAfragmat = sys.argv[2] # args.DNAfragmat
gene_data = sys.argv[3] # args.gene_info
vcf = sys.argv[4] #args.vcf
outputname = sys.argv[5] #args.outputfolder
isoforms = sys.argv[6] # args.isoforms

time_start = time.time()

#if vcf == None:
#    raise ValueError("Error: vcf required")

# if RNAfragmat == None and DNAfragmat == None:
#     sys.exit(
#         "Error: at least one of RNAfragmat or DNAfragmat should be given. You can use Chair tool to generate them from SAM files."
#     )

RNAfragmats = [RNAfragmat]
DNAfragmats = [DNAfragmat]

if not os.path.exists(outputname):
    os.system (f"mkdir -p {outputname}")
pair_thresh = 0.7


def make_golden_from_true2(
    filename: str,
) -> tuple[dict[int, dict[int, str]], list[str], list[str]]:
    f = open(filename, "r")
    a = f.readlines()
    f.close()
    i = 0
    while a[i][0] == "#":
        i += 1
    c = [x.split()[9][:3] for x in a[i:]]
    vChroms = [x.split()[0] for x in a[i:]]
    vPositions = [x.split()[1] for x in a[i:]]
    d = {0: {}, 1: {}}
    for j in range(len(c)):
        if c[j][1] == "|":
            d[0][j] = int(c[j][0])
            d[1][j] = int(c[j][2])
        else:
            d[0][j] = "."
            d[1][j] = "."
    return d, vChroms, vPositions
print("Loading VCF file for scoring")
V, vcfChroms, vcfPositions = make_golden_from_true2(vcf)

if RNAfragmat == None or gene_data == None:
    print("Missing RNAfragmat, running phasing without DASE")
    if DNAfragmat == None:  # Only RNA data provided without gene info
        D = datagen.make_data_from_fragmat(RNAfragmats, vcf, 0.02)
        G = basic_class.easy_graph(D)
        GX = alg.RNA_phase(
            0.001, pair_thresh, 0.02, G.read_dict, G.comp_mins, G.components
        )
        output.make_solution(
            GX, G, outputname, "HapTreeX_noDASE_output.txt"
        )  # Running regular Haptree on 2-reads in the genome
    else:
        if RNAfragmat != None:
            DNAfragmats = [DNAfragmat, RNAfragmat]
        D = datagen.make_data_from_fragmat(DNAfragmats, vcf, 0.02)
        time1 = time.time()
        print(("Pt 1 took ", time1 - time_start))
        G = basic_class.easy_graph(D)
        time2 = time.time()
        print(("Pt 2 took ", time2 - time1))
        GX = alg.RNA_phase(
            0.001, pair_thresh, 0.02, G.read_dict, G.comp_mins, G.components
        )
        time3 = time.time()
        print(("Pt 3 took ", time3 - time2))
        output.make_solution(
            GX, G, outputname, "HapTreeX_noDASE_output.txt"
        )  # Running regular Haptree on 2-reads in the genome
        time4 = time.time()
        print(("Pt 4 took ", time4 - time3))
        print(("Total took ", time4 - time_start))
else:
    RD = datagen.make_RNA_data_from_fragmat(gene_data, RNAfragmats, vcf, 0.02, isoforms)
    time1 = time.time()
    print(("Pt 1 took ", time1 - time_start))
    stats.stats(RD, 2, 0.2, 0.6, 0, 2, 0.001, 0.2)
    times = time.time()
    print(("stats took ", times - time1))
    if DNAfragmat == None:
        D = rna_class.make_data_from_RNA_data(RD)
    else:
        D = datagen.make_data_from_fragmat(DNAfragmats, vcf, 0.02, RD.long_read_list)
    G = basic_class.easy_graph(D)
    jG = joint_class.joint_graph(RD, G)
    time2 = time.time()
    print(("Pt 2 took ", time2 - times))
    jGX = alg.RNA_phase(
        0.001, pair_thresh, 0.02, jG.read_dict, jG.comp_mins, jG.components
    )
    time3 = time.time()
    print(("Pt 3 took ", time3 - time2))
    output.make_solution_simple(
        jGX,
        outputname,
        "HapTreeX_withDASE_output.txt",
        vcfChroms,
        vcfPositions,
    )  # Phasing using both 1-reads and 2-reads from DNAfragmat + RNAfragmat
    time4 = time.time()
    print(("Pt 4 took ", time4 - time3))
    print(("Total took ", time4 - time_start))
