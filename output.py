import subprocess

from typing import dict, list


def make_solution(X, G, outputname, name):
    # writes phasing solution to outputname
    t = "\t"
    f = open(outputname + "/" + name, "w")
    for start in G.comp_mins:
        s_pos = G.nodes[start].position
        e_pos = G.nodes[G.components[start][-1]].position
        f.write(
            "BLOCK"
            + " Start: "
            + str(start + 1)
            + " Len: "
            + str(len(G.components[start]))
            + " Phased: "
            + str(len(G.components[start]))
            + " Span: "
            + str(e_pos - s_pos)
            + " MEC: "
            + str(comp_MEC(start, X, G))
            + " Reads: "
            + str(len(G.comp_reads[start]))
            + " \n"
        )

        for i in G.components[start]:
            s = ""
            s = (
                s
                + str(i + 1)
                + t
                + str(X[start][0][i])
                + t
                + str(X[start][1][i])
                + t
                + str(G.nodes[i].chrom)
                + t
                + str(G.nodes[i].position)
                + t
            )
            s += "\n"
            f.write(s)
        f.write("*****\n")
    f.close()


def make_solution_RNA(X, RD, outputname, name):
    # writes phasing solution to outputname
    t = "\t"
    f = open(outputname + "/" + name, "w")
    for start in X:
        end = max(X[start][0].keys())
        s_pos = RD.nodes[start].position
        e_pos = RD.nodes[end].position
        f.write(
            "BLOCK"
            + " Start: "
            + str(start + 1)
            + " Len: "
            + str(len(X[start][0]))
            + " Phased: "
            + str(len(X[start][0]))
            + " Span: "
            + str(e_pos - s_pos)
            + " MEC: "
            + "##"
            + " Reads: "
            + "##"
            + " \n"
        )

        for i in X[start][0]:
            s = ""
            s = (
                s
                + str(i + 1)
                + t
                + str(X[start][0][i])
                + t
                + str(X[start][1][i])
                + t
                + str(RD.nodes[i].chrom)
                + t
                + str(RD.nodes[i].position)
                + t
            )
            s += "\n"
            f.write(s)
        f.write("*****\n")
    f.close()


def make_solution_simple(
    X: dict[int, dict[int, dict[int, int]]],
    outputname: str,
    name: str,
    vChroms: list[str],
    vPositions: list[str],
) -> None:
    # writes phasing solution to outputname
    t = "\t"
    f = open(outputname + "/" + name, "w")
    for start in X:
        f.write(
            "BLOCK"
            + " Start: "
            + str(start + 1)
            + " Len: "
            + str(len(X[start][0]))
            + " Phased: "
            + str(len(X[start][0]))
            + " Span: "
            + "##"
            + " MEC: "
            + "##"
            + " Reads: "
            + "##"
            + " \n"
        )

        for i in X[start][0]:
            s = ""
            s = (
                s
                + str(i + 1)
                + t
                + str(X[start][0][i])
                + t
                + str(X[start][1][i])
                + t
                + str(vChroms[i])
                + t
                + str(vPositions[i])
                + t
            )
            s += "\n"
            f.write(s)
        f.write("*****\n")
    f.close()


def comp_MEC(m, X, G):
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


def make_golden_from_true(filename, G):
    f = open(filename, "r")
    a = f.readlines()
    f.close()
    c = [x.split("\t")[10][:3] for x in a]
    d = {m: {0: {}, 1: {}} for m in G.comp_mins}
    for m in G.comp_mins:
        for j in G.components[m]:
            if c[j][1] == "|":
                d[m][0][j] = int(c[j][0])
                d[m][1][j] = int(c[j][2])
            else:
                d[m][0][j] = "."
                d[m][1][j] = "."
    return d
