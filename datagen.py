from gene import Gene
from read import Read
from rna import RNAData
from graph import Data, edges_from_readlist
from common import QUALITY_CUTOFF


###########################################################################
# GENES

def determine_genes_gtf(gene_data: str, chroms: set[str]) -> dict[int, Gene]:
    print("Loading and formatting genes")
    f = open(gene_data, "r")
    a = list(f.readlines())
    f.close()
    i = 0
    while a[i][0] == "#":
        i += 1
    c = a[i:]
    b = dict[int, list[str]]()
    for j in range(len(a[i:])):
        b[j] = c[j][:-2].split("\t")
    trans_list = list[int]()
    trans_dict = dict[int, dict[str, list[int]]]()
    for j in range(len(b)):
        y = b[j]
        if y[2] == "transcript":
            trans_list.append(j)
            trans_dict[j] = {"exon_starts": [], "exon_lengths": []}
    for k in range(len(trans_list) - 1):
        start = trans_list[k]
        end = trans_list[k + 1]
        for j in range(start, end):
            y = b[j]
            if y[2] == "exon":
                exon_start = int(y[3])
                exon_end = int(y[4])
                exon_length = exon_end - exon_start + 1
                trans_dict[start]["exon_starts"].append(exon_start)
                trans_dict[start]["exon_lengths"].append(exon_length)

    start = trans_list[-1]
    end = len(b)
    for j in range(start, end):
        y = b[j]
        if y[2] == "exon":
            exon_start = int(y[3])
            exon_end = int(y[4])
            exon_length = exon_end - exon_start + 1
            trans_dict[start]["exon_starts"].append(exon_start)
            trans_dict[start]["exon_lengths"].append(exon_length)

    transcripts = dict[int, Gene]()
    for j in trans_dict:
        y = b[j]
        chrom = y[0]
        if chrom in chroms:
            bp_start = int(y[3])
            bp_end = int(y[4])
            sign = y[6]
            z = y[8]
            w = z.split("; ")
            v = [xx.split(" ") for xx in w]
            gene_id = v[0][1][1:-1]
            transcript_id = v[1][1][1:-1]
            gene_type = v[2][1][1:-1]
            exon_starts = trans_dict[j]["exon_starts"]
            exon_count = len(trans_dict[j]["exon_lengths"])
            exon_lengths = trans_dict[j]["exon_lengths"]
            transcripts[j] = Gene(
                transcript_id,
                chrom,
                bp_start,
                bp_end,
                sign,
                exon_count,
                exon_lengths,
                exon_starts,
                gene_id,
                gene_type,
            )
            transcripts[j].index = j

    return transcripts


def build_isodict(isoforms: str) -> dict[str, list[float, str]]:
    isodict = dict[str, list[float, str]]()
    f = open(isoforms)
    temp = [x.split("\t") for x in f.readlines()]
    for i in range(1, len(temp)):
        isodict[temp[i][0]] = [float(temp[i][9]), temp[i][3]]
    return isodict


def filter_transcripts(genes, isodict: dict[str, list[float, str]]) -> dict[int, Gene]:
    gene_cov = dict[str, int]()
    for i in isodict:
        gene_cov[isodict[i][1]] = 0
    for i in isodict:
        gene_cov[isodict[i][1]] += isodict[i][0]

    new_genes = dict[int, Gene]()
    types = set[str]()
    for j in genes:
        g = genes[j]
        tid = g.transcript_id
        types.add(genes[j].gene_type)
        if tid in isodict:
            val = isodict[tid][0]
            if val > 0:
                if val / gene_cov[isodict[tid][1]] > 0.01:
                    new_genes[j] = g
    return new_genes


###########################################################################
# READS


def make_read_of_frag(frag0: str) -> tuple[list[tuple[int, int]], int]:
    frag = frag0.split(" ")[2:]
    # takes line from fragment matrix and makes tuple formatted read
    read = list[tuple[int, int]]()
    qual = frag[-1][0]
    frag = frag[:-1]
    if len(frag) % 2 == 1:
        raise ValueError(f"fragment file error: {frag}")

    for i in range(0, len(frag), 2):
        key = int(frag[i])
        for char in frag[i + 1]:
            read.append((key - 1, int(char)))
            key += 1
    if len(read) == 1:
        qual = ord(qual)
    else:
        qual = 1000
    return read, qual


def make_readlist_from_fragmat(
    fragmats: list[str],
    skip_single: bool = False
) -> dict[int, Read]:
    # translates fragment matrix into a list of reads
    print("Loading and formatting fragments")
    a = list[str]()
    for fragmat in fragmats:
        f = open(fragmat, "r")
        A = list(f.readlines())
        f.close()
        a = a + A

    F = len(a)
    read_list_list = list[list[tuple[int, int]]]()
    i = 0
    for r in a:
        i += 1
        read, qual = make_read_of_frag(r)
        if qual >= QUALITY_CUTOFF:
            if not skip_single or len(read) > 1:
                read_list_list.append(read)

    print(f"{len(read_list_list)} reads of sufficient quality")
    read_list = dict[int, Read]()
    read_counter = dict[list[tuple[int, int]], int]()
    i = 0
    for tup_read in read_list_list:
        i += 1
        if tup_read in read_counter:
            read_counter[tup_read] += 1
        else:
            read_counter[tup_read] = 1
    i = 0
    print(f"{len(read_counter)} distinct reads")
    for tup in read_counter:
        read = dict[list[tuple[int, int]], int]()
        for k, v in tup:
            read[k] = v
        read_list[i] = Read(read, read_counter[tup], i)
        i += 1
    return read_list


###########################################################################
# make RNA_data and DNA_data objects


def positions_names_states(
    vcf: str,
) -> tuple[dict[int, int], dict[int, str], dict[int, str], dict[int, int], int]:
    # reading data from VCF file and formatting as lists
    f = open(vcf, "r")
    a = f.readlines()
    i = 0
    while a[i][0] == "#":
        i += 1
    f.close()
    b = [x.split("\t") for x in a[i:]]
    positions = dict[int, int]()
    states = dict[int, int]()
    names = dict[int, str]()
    chroms = dict[int, str]()
    for i in range(len(b)):
        line = b[i]
        positions[i] = int(line[1])
        names[i] = line[2]  ##change when fixed files
        chroms[i] = line[0]
        states[i] = 1  # sum(line[9]) ??
    k = 2
    return (states, names, chroms, positions, k)


def make_RNA_data_from_fragmat(
    gene_data: str,
    fragmats: list[str],
    vcf: str,
    error: float,
    isoforms: str
) -> RNAData:
    ##RNA DATA
    read_list = make_readlist_from_fragmat(fragmats)
    print(f"Loading VCF file {vcf}")
    S, names, chroms, positions, k = positions_names_states(vcf)
    chrom_set = set(chroms.values())
    n = len(S)
    print("Preparing data for ReadGraph")
    genes = determine_genes_gtf(gene_data, chrom_set)

    isodict = None
    filtered_genes = genes
    if isoforms != "":
        print("Building IsoDict")
        isodict = build_isodict(isoforms)
        filtered_genes = filter_transcripts(genes, isodict)

    return RNAData(
        S, genes, filtered_genes, error, read_list, positions, names, chroms, isodict
    )


def make_data_from_fragmat(
    fragmat: list[str],
    vcf: string,
    error: float,
    RNA_readlist: dict[int, Read]
) -> Data:
    ##regular DNA fragmat
    read_list = make_readlist_from_fragmat(fragmat, skip_single=True)
    if len(RNA_readlist) > 0:
        max_key = max(read_list.keys())
        for zz in RNA_readlist:
            read_list[zz + max_key] = RNA_readlist[zz]
    for r in read_list.values():
        r.special_key = r.keys[1]
        r.rates = {0: 0.5, 1: 0.5}

    print("Loading VCF file")
    S, names, chroms, positions, k = positions_names_states(vcf)
    n = len(S)
    print(f"{n} SNPs in VCF file")

    D = edges_from_readlist(read_list)
    return Data(D, S, k, error, read_list, positions, names, chroms)
