from basic_class import Graph
from rna_class import RNAData


def comp_MEC(m: int, X: dict[int, dict[int, dict[int, int]]], G: Graph) -> int:
    # computes MEC score of component containing node m
    total = 0
    for read in G.comp_reads[m]:
        counts = [0] * G.k
        for key in read.keys:
            for i in range(G.k):
                if not X[m][i][key] == read.read[key]:
                    counts[i] += 1
        total += min(counts)
    return total


def make_solution(
    X: dict[int, dict[int, dict[int, int]]],
    G: Graph,
    outputname: str,
    name: str
):
    # writes phasing solution to outputname
    with open(f"{outputname}/{name}", "w") as f:
        for start in G.comp_mins:
            s_pos = G.nodes[start].position
            e_pos = G.nodes[G.components[start][-1]].position
            f.write(
                f"BLOCK Start: {start + 1} Len: {len(G.components[start])} "
                + f"Phased: {len(G.components[start])} Span: {e_pos - s_pos} "
                + f"MEC: {comp_MEC(start, X, G)} Reads: {len(G.comp_reads[start])}\n"
            )
            for i in G.components[start]:
                f.write(
                    f"{i+1}\t{X[start][0][i]}\t{X[start][1][i]}\t{G.nodes[i].chrom}\t"
                    + f"{G.nodes[i].position}\t\n"
                )
            f.write("*****\n")


def make_solution_RNA(
    X: dict[int, dict[int, dict[int, int]]],
    RD: RNAData,
    outputname: str,
    name: str
):
    # writes phasing solution to outputname
    with open(f"{outputname}/{name}", "w") as f:
        for start in X:
            end = max(X[start][0].keys())
            s_pos = RD.nodes[start].position
            e_pos = RD.nodes[end].position

            f.write(
                f"BLOCK Start: {start + 1} Len: {len(X[start][0])} "
                + f"Phased: {len(X[start][0])} Span: {e_pos - s_pos} "
                + f"MEC: ## Reads: ##\n"
            )
            for i in X[start][0]:
                f.write(
                    f"{i+1}\t{X[start][0][i]}\t{X[start][1][i]}\t" +
                    f"{RD.G.nodes[i].chrom}\t{RD.G.nodes[i].position}\t\n"
                )
            f.write("*****\n")


def make_solution_simple(
    X: dict[int, dict[int, dict[int, int]]],
    outputname: str,
    name: str,
    vChroms: list[str],
    vPositions: list[str],
):
    # writes phasing solution to outputname
    with open(f"{outputname}/{name}", "w") as f:
        for start in X:
            f.write(
                f"BLOCK Start: {start + 1} Len: {len(X[start][0])} "
                + f"Phased: {len(X[start][0])} Span: ## " +
                + f"MEC: ## Reads: ##\n"
            )
            for i in X[start][0]:
                f.write(
                    f"{i+1}\t{X[start][0][i]}\t{X[start][1][i]}\t{vChroms[i]}\t"
                    + f"{vPositions[i]}\t\n"
                )
            f.write("*****\n")

