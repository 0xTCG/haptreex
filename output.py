from read import SNP, Read
from datagen import VCF
from alg import Phase
from graph import Graph
# from rna import RNAData
from typing import Tuple, Dict, List, Set


def mec_score(ploidy: int, reads: List[Read], phase: Phase) -> int:
    """
    Computes a phase MEC score of a connected component.
    """
    return sum(
        min(
            sum(1 for snp in read.snps if phase.haplotypes[hap][snp] != read.snps[snp])
            for hap in range(ploidy)
        )
        for read in reads
    )


def make_solution(V: VCF, G: Graph, phases: Dict[int, Phase], path: str):
    # polyploid-ready
    with open(path, "w") as f:
        for root, comp in sorted(G.components.items()):
            span = V.snps[comp.nodes[-1]].pos - V.snps[comp.nodes[0]].pos
            mec = mec_score(G.ploidy, comp.reads, phases[root])
            reads = sum(len(x) for x in comp.reads)
            span=mec=reads="##"
            f.write(
                f"BLOCK "
                f"Start: {root + 1} "
                f"Len: {len(comp.nodes)} "
                f"Phased: {len(comp.nodes)} "
                f"Span: {span} "
                f"MEC: {mec} "
                f"Reads: {reads}\n"
            )
            haps = sorted(phases[root].haplotypes, key=lambda d: d[comp.nodes[0]])
            for snp in comp.nodes:
                f.write(f"{snp + 1}\t")
                for hap in range(G.ploidy):
                    f.write(f"{haps[hap][snp]}\t")
                f.write(f"{V.snps[snp].chr}\t{V.snps[snp].pos + 1}\t\n")
            f.write("*****\n")


# def make_solution_RNA(
#     X: Dict[int, Tuple[Dict[int, int], Dict[int, int]]],
#     RD: RNAData,
#     outputname: str
# ):
#     # writes phasing solution to outputname
#     with open(outputname, "w") as f:
#         for start in sorted(X):
#             end = max(X[start][0].keys())
#             s_pos = RD.data.nodes[start].pos
#             e_pos = RD.data.nodes[end].pos
#             f.write(
#                 f"BLOCK Start: {start + 1} Len: {len(X[start][0])} "
#                 + f"Phased: {len(X[start][0])} Span: {e_pos - s_pos} "
#                 + f"MEC: ## Reads: ##\n"
#             )
#             for i in sorted(X[start][0]):
#                 f.write(
#                     f"{i+1}\t{X[start][0][i]}\t{X[start][1][i]}\t" +
#                     f"{RD.G.data.nodes[i].chrom}\t{RD.G.data.nodes[i].pos}\t\n"
#                 )
#             f.write("*****\n")


# def make_solution_simple(
#     X: Dict[int, Tuple[Dict[int, int], Dict[int, int]]],
#     outputname: str,
#     vChroms: List[str],
#     vPositions: List[str]
# ):
#     # writes phasing solution to outputname
#     with open(outputname, "w") as f:
#         for start in sorted(X):
#             f.write(
#                 f"BLOCK Start: {start + 1} Len: {len(X[start][0])} "
#                 + f"Phased: {len(X[start][0])} Span: ## " +
#                 + f"MEC: ## Reads: ##\n"
#             )
#             for i in sorted(X[start][0]):
#                 f.write(
#                     f"{i+1}\t{X[start][0][i]}\t{X[start][1][i]}\t{vChroms[i]}\t"
#                     + f"{vPositions[i]}\t\n"
#                 )
#             f.write("*****\n")
