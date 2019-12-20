from read import Read
from files import VCF
from alg import Phase
from graph import Graph
from typing import Dict, List


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
            # span = mec = reads = "##"  # For debugging purposes
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
