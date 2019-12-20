import files
import graph
import output
import sys
import alg
import random
# import chair
from getopt import gnu_getopt as getopt, GetoptError
from typing import Tuple

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


vcf_path, dna, rna, gtf, out, isoforms = parse_args()

PHASE_THRESHOLD = 0.001
PAIR_THRESHOLD = 0.7
PHASE_ERROR = 0.02

DASE_ONLY = True

print(f"Loading VCF file {vcf_path}...")
vcf = files.parse_vcf(vcf_path)
print(f"{len(vcf.snps)} SNPs in VCF file")

if dna == "" and rna == "":
    print("Need either RNA or DNA data to phase")
elif gtf == "" or rna == "":
    print("Running DNA/RNA phasing without DASE (no gene data provided)")
    dnas = ([dna] if dna else []) + ([rna] if rna else [])
    g = files.load_dna_data(vcf, dnas)
else:
    r = files.load_rna_data(vcf, gtf, [rna], isoforms)
    if not dna:
        g = graph.Graph(r.multi_reads if not DASE_ONLY else [], r.ploidy)
    else:
        g = files.load_dna_data(vcf, [dna], r.multi_reads)
    g.integrate_rna(r)
phases = alg.phase(g, PHASE_THRESHOLD, PAIR_THRESHOLD, PHASE_ERROR)
output.make_solution(vcf, g, phases, out)
