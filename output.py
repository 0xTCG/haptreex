from graph import Graph
from rna import RNAData
from typing import Tuple, Dict, List, Set


def mec_score(ploidy: int, reads: List[Read], phase: Phase) -> int:
    """
    Computes a phase MEC score of a connected component.
    """
    total = 0
    for read in reads:
        counts = [0] * ploidy
        for snp in read.snps:
            for hap in range(ploidy):
                if not phase.haplotypes[hap][snp] == read.snps[snp]:
                    counts[hap] += 1
        total += min(counts)
    return total


def make_solution(G: Graph, phases: Dict[SNP, Phase], path: str):
    # polyploid-ready
    with open(path, "w") as f:
        for root in G.component_roots:
            s_pos = G.data.nodes[root].snp.pos
            e_pos = G.data.nodes[G.components[root][-1]].snp.pos
            mec = mec_score(G.ploidy, G.component_reads[root], phases[root])
            reads = sum(len(x) for x in G.component_reads[root])
            f.write(
                f"BLOCK "
                f"Start: {root} "
                f"Len: {len(G.components[root])} "
                f"Phased: {len(G.components[root])} "
                f"Span: {e_pos - s_pos} "
                f"MEC: {mec} "
                f"Reads: {reads}\n"
            )
            for snp in G.components[root]:
                f.write(f"{snp}\t")
                for hap in range(G.ploidy):
                    f.write(f"{phases[root].haplotypes[hap][snp]}\t")
                f.write(f"{snp.chr}\t{snp.pos}\t\n")
            f.write("*****\n")


def make_solution_RNA(
    X: Dict[int, Tuple[Dict[int, int], Dict[int, int]]],
    RD: RNAData,
    outputname: str
):
    # writes phasing solution to outputname
    with open(outputname, "w") as f:
        for start in sorted(X):
            end = max(X[start][0].keys())
            s_pos = RD.data.nodes[start].pos
            e_pos = RD.data.nodes[end].pos
            f.write(
                f"BLOCK Start: {start + 1} Len: {len(X[start][0])} "
                + f"Phased: {len(X[start][0])} Span: {e_pos - s_pos} "
                + f"MEC: ## Reads: ##\n"
            )
            for i in sorted(X[start][0]):
                f.write(
                    f"{i+1}\t{X[start][0][i]}\t{X[start][1][i]}\t" +
                    f"{RD.G.data.nodes[i].chrom}\t{RD.G.data.nodes[i].pos}\t\n"
                )
            f.write("*****\n")


def make_solution_simple(
    X: Dict[int, Tuple[Dict[int, int], Dict[int, int]]],
    outputname: str,
    vChroms: List[str],
    vPositions: List[str]
):
    # writes phasing solution to outputname
    with open(outputname, "w") as f:
        for start in sorted(X):
            f.write(
                f"BLOCK Start: {start + 1} Len: {len(X[start][0])} "
                + f"Phased: {len(X[start][0])} Span: ## " +
                + f"MEC: ## Reads: ##\n"
            )
            for i in sorted(X[start][0]):
                f.write(
                    f"{i+1}\t{X[start][0][i]}\t{X[start][1][i]}\t{vChroms[i]}\t"
                    + f"{vPositions[i]}\t\n"
                )
            f.write("*****\n")

