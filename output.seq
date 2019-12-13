from graph import Graph
from rna import RNAData


def comp_MEC(m: int, X: dict[int, tuple[dict[int, int], dict[int, int]]], G: Graph) -> int:
    # computes MEC score of component containing node m
    # polyploid-ready
    total = 0
    assert G.data.k == 2
    for read in G.comp_reads[m]:
        counts = [0] * G.data.k
        for key in read.keys:
            for i in range(G.data.k):
                if not X[m][i][key] == read.read[key]: 
                    counts[i] += 1
        total += min(counts)
    return total


def make_solution(
    X: dict[int, tuple[dict[int, int], dict[int, int]]],
    G: Graph,
    outputname: str
):
    # writes phasing solution to outputname
    # polyploid-ready
    with open(outputname, "w") as f:
        for start in sorted(G.comp_mins):
            s_pos = G.data.nodes[start].pos
            e_pos = G.data.nodes[G.components[start][-1]].pos
            reads = sum(x.count for x in G.comp_reads[start])
            f.write(
                f"BLOCK Start: {start + 1} Len: {len(G.components[start])} "
                + f"Phased: {len(G.components[start])} Span: {e_pos - s_pos} "
                + f"MEC: {comp_MEC(start, X, G)} Reads: {reads}\n"
            )
            for i in G.components[start]:
                f.write(f"{i+1}\t")
                for p in range(G.data.k):
                    f.write(f"{X[start][p][i]}\t")
                f.write(f"{G.data.nodes[i].chrom}\t{G.data.nodes[i].pos}\t\n")
            f.write("*****\n")


def make_solution_RNA(
    X: dict[int, tuple[dict[int, int], dict[int, int]]],
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
    X: dict[int, tuple[dict[int, int], dict[int, int]]],
    outputname: str,
    vChroms: list[str],
    vPositions: list[str]
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

