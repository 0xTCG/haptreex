from read import SNP, Read
from graph import Graph
from common import QUALITY_CUTOFF
from typing import Tuple, Dict, List, Set, NamedTuple, Optional, Iterator, Any
from dataclasses import dataclass
from RNA import RNAData
from gene import Gene
import pysam
import bisect
from pprint import pprint


@dataclass
class VCF:
    snps: List[SNP]
    line_to_snp: Dict[int, int] # 1-based index
    chromosomes: Set[str]


def parse_vcf(
    vcf_path: str,
    sample: Optional[str] = None
) -> VCF:
    """
    Parse a VCF file and return a dictionary of sorted SNPs:
        d := {chromosome_name: sorted([snp_1, snp_2])}
    and the corresponding line index
        {line_in_vcf: (chr, index in d[chr])}.
    Each snp_i is a SNP that is heterozygous in the sample (i.e. |set(GT(snp_i))| > 1).
    """

    snps: List[SNP] = []
    line_to_snp: Dict[int, int] = {}
    chromosomes: Set[str] = set()
    with open(vcf_path) as vcf:
        # samples = list(vcf.header.samples)
        # if not samples:
        #     raise ValueError('No samples present in the VCF file')
        # if sample and sample not in samples:
        #     raise ValueError(f'Sample {sample} not found in the VCF')
        # sample = sample if sample else samples[0]
        print(f'Using sample {sample} in {vcf_path}...')
        seen_snps = 0
        for line in vcf:
            if line[0] == '#':
                continue
            seen_snps += 1
            chr, pos, id, ref, alt, _, _, _, fmt_, sample = line.split('\t')
            if len(ref) != 1:  # Ignore indels
                continue
            # We only deal with SNPs here for now
            fmt = dict(zip(fmt_.split(':'), sample.split(':')))
            gt = fmt['GT'].replace('|', '/').split('/')
            alleles = [ref] + alt.split(',')
            # Get only alleles that are specified in GT field
            alleles = [a for i, a in enumerate(alleles) if str(i) in gt and len(a) == 1]
            if len(alleles) > 1:  # Ignore homozygous SNPs
                snp = SNP(len(snps), chr, int(pos) - 1, id, alleles)
                if snps and snp < snps[-1]:  # TODO
                    raise ValueError(f'VCF is not sorted (SNP {snp})')
                line_to_snp[seen_snps] = len(snp)
                snps.append(snp)
                chromosomes.add(chr)

    return VCF(snps, line_to_snp, chromosomes)


def parse_gtf(gtf_path: str, chroms: Set[str]) -> Iterator[Gene]:
    b = [a[:-2].split("\t") for a in open(gtf_path, "r") if a[0] != '#']

    genes: List[int] = []
    exons: List[List[Tuple[int, int]]] = []  # Gene -> List of exon (start, len)
    for gene, line in enumerate(b):
        if line[2] == "transcript":
            genes.append(gene)
            exons.append([])
    for gene in range(len(genes)):
        end = genes[gene + 1] if gene < len(genes) - 1 else len(b)
        for l in range(genes[gene], end):
            if b[l][2] == "exon":
                exon_start = int(b[l][3])
                exon_end = int(b[l][4])
                exons[gene].append((exon_start, exon_end - exon_start + 1))
    for gene in genes:
        chrom = b[gene][0]
        if chrom not in chroms:
            continue
        v = [xx.split(" ") for xx in b[gene][8].split("; ")]
        yield Gene(
            id=gene,
            name=v[1][1][1:-1],
            chr=chrom,
            interval=(int(b[gene][3]), int(b[gene][4])),
            sign=b[gene][6],
            exons=exons[gene],
        )


def parse_read(
    lines,
    snps: VCF,
    threshold: int,
    ignore_conflicts: bool = True
) -> Iterator[Tuple[str, List[Tuple[SNP, int, str]]]]:
    """
    If reads are valid and pass the threshold filter, yields the
        (read_name, [all_1, all_2, ...])
    where
        all_i := (snp, allele, quality).
    Example:
        ('read1', [(SNP("chr1", 12), 'A', 'E'), ...])
    """

    cov: Dict[SNP, Dict[int, List[str]]] = {}  # SNP: {allele: [qual1, qual2, ...]}
    name = lines[0].query_name
    counts = [0] * len(lines)
    for line_i, sam in enumerate(lines):
        read, ref = 0, 0
        cand = snps[sam.reference_name]
        for op, sz in sam.cigartuples:
            if op in [0, 7, 8]:  # 'M=X':
                s = SNP(sam.reference_name, sam.pos + ref, "", [])
                x = bisect.bisect_left(cand, s)
                for i in range(x, len(cand)):
                    snp = cand[i]
                    if snp.pos >= sam.pos + ref + sz:
                        break
                    t = snp.pos - sam.pos - ref + read
                    if sam.seq[t] in snp.alleles:
                        allele = snp.alleles.index(sam.seq[t])
                        qual = sam.qual[t] if sam.qual != '*' else '.'
                        cov.setdefault(snp, {}).setdefault(allele, []).append(qual)
                        counts[line_i] += 1
                read += sz
                ref += sz
            elif op in [1, 4]:  # 'IS':
                read += sz
            elif op in [2, 3, 5, 6]:  # 'DNHP':
                ref += sz
    # Adding an MP suffix like extractHairs to specify that matepairs are merged
    if len(counts) > 1 and counts[0] > 0 and counts[1] > 0:
        name += "_MP"
    # Filter out SNPs that harbour mate-pair allele conflicts
    for snp in list(cov):
        if ignore_conflicts and len(cov[snp]) > 1:
            del cov[snp]
    if len(cov) >= threshold:
        yield name, [(snp, a, max(q)) for snp, A in cov.items() for a, q in A.items()]


def parse_bam(
    sam_path: str,
    snps: VCF,
    threshold: int = 1,
    no_chimeric: bool = True,
    no_duplicates: bool = True,
    no_conflicts: bool = True
) -> Iterator[Tuple[str, List[Tuple[SNP, int, str]]]]:
    """
    Reads a sorted SAM/BAM.
    """

    seen: Dict[str, Any] = {}  #S Dict[str, SAMRecord]
    seen_chrs: Set[str] = set()
    with pysam.AlignedSegment(sam_path) as sam:
        for line in sam:
            if line.reference_name not in seen_chrs:
                print(f'Parsing {line.reference_name}, {len(seen)} cached so far...')
            seen_chrs.add(line.reference_name)
            name = sam.query_name
            if len(name) > 2 and name[-2] in '#/':
                name = name[:-2]

            if (
                line.is_supplementary
                or (no_duplicates and line.is_duplicate)
                or (no_chimeric and line.reference_name != line.next_reference_name)
                or line.is_unmapped
                or not line.cigartuples
                or line.reference_name not in snps.chromosomes
            ):
                if (
                    not line.is_supplementary
                    and not line.mate_is_unmapped
                    and name in seen
                ):
                    yield from parse_read([seen[name]], snps, threshold, no_conflicts)
                    del seen[name]
                continue
            elif name in seen:
                yield from parse_read([seen[name], line], snps, threshold, no_conflicts)
                del seen[name]
            elif (
                not line.mate_is_unmapped
                and line.reference_name == line.next_reference_name
                and line.mpos < line.pos
            ):
                yield from parse_read([line], snps, threshold, no_conflicts)
            elif (
                not line.mate_is_unmapped
                and line.reference_name != line.next_reference_name
                and line.next_reference_name in seen_chrs
            ):
                yield from parse_read([line], snps, threshold, no_conflicts)
            else:
                seen[name] = line
    for line in seen:
        yield from parse_read([line], snps, threshold, no_conflicts)


def parse_fragmat(
    fragmat: str,
    vcf: VCF,
    skip_single: bool
) -> Iterator[Tuple[str, List[Tuple[int, int, str]]]]:
    """
    Yields reads in the format
        (read_name, [(snp_id, allele_id, qual), ...])
    """
    print(f"Loading and formatting fragments...")
    with open(fragmat, "r") as f:
        for r in f:
            frag = r.split()
            name, qual, frag = frag[1], frag[-1], frag[2:-1]
            if len(frag) % 2 == 1:
                raise ValueError(f"fragment file error: {frag}")
            if len(qual) == 1 and ord(qual) < QUALITY_CUTOFF:  # TODO: fix this
                continue

            alleles: List[Tuple[int, int, str]] = []
            for i in range(0, len(frag), 2):
                idx = int(frag[i])
                for allele in frag[i + 1]:
                    if idx not in vcf.line_to_snp:
                        print(f'Warning: Invalid SNP index {idx} for read {name}')
                        continue
                    snp = vcf.snps[vcf.line_to_snp[idx]]
                    if not 0 <= int(allele) < len(snp.alleles):
                        raise ValueError(f'Invalid allele {allele} for SNP {snp}')
                    alleles.append((idx, int(allele), qual[len(alleles)]))
                    idx += 1
            yield name, alleles


def parse_phases(
    vcf: VCF,
    paths: List[str],
    skip_single: bool = True
) -> Iterator[Read]:
    reads = []
    for path in paths:
        print(f'Parsing {path}...')
        for _, alleles in parse_fragmat(path, vcf, skip_single):
            if not (skip_single and len(alleles) <= 1):
                reads.append(alleles)
    print(f"{len(reads)} reads of sufficient quality")

    read_counter = {}  #S : Dict[List[Tuple[int, int]], int] = {}
    for r in reads:
        key = tuple((s, a) for s, a, _ in r)
        if key not in read_counter:
            read_counter[key] = 1
        else:
            read_counter[key] += 1
    print(f"{len(read_counter)} distinct reads")

    for i, (tup, cnt) in enumerate(read_counter.items()):
        yield Read({snp_id: snp_a for snp_id, snp_a in tup}, cnt, i)


###########################################################################
# make RNA_data and DNA_data objects


def load_rna_data(
    vcf_path: str,
    gtf_path: str,
    paths: List[str],
    error: float,
    isoforms_path: str
) -> RNAData:
    print(f"Loading VCF file {vcf_path}...")
    vcf = parse_vcf(vcf_path)
    print(f"{len(vcf.snps)} SNPs in VCF file")

    print(f"Loading GTF {gtf_path}...")
    genes = list(parse_gtf(gtf_path, vcf.chromosomes))
    print(f"{len(genes)} genes in GTF file")

    # graph = Graph(reads, ploidy=2, error=error)
    # return vcf, graph
    isodict: Dict[str, Tuple[float, str]] = {}
    filtered_genes = genes
    if isoforms != "":
        print("Building IsoDict")
        isodict = build_isodict(isoforms)
        filtered_genes = filter_transcripts(genes, isodict)

    return RNAData(vcf, filtered_genes, reads, error)


def load_dna_data(
    vcf_path: str,
    paths: List[str],
    error: float,
    RNA_readlist: Dict[int, Read] = None
) -> Tuple[VCF, Graph]:
    print(f"Loading VCF file {vcf_path}...")
    vcf = parse_vcf(vcf_path)
    print(f"{len(vcf.snps)} SNPs in VCF file")

    reads = list(parse_phases(vcf, paths, skip_single=True))
    # if len(RNA_readlist) > 0:
    #     max_key = max(reads.keys())
    #     for zz in RNA_readlist:
    #         reads[zz + max_key] = RNA_readlist[zz]
    for r in reads:
        r.special_snp = sorted(r.snps)[1]
        r.rates = [0.5, 0.5]

    graph = Graph(reads, ploidy=2, error=error)
    return vcf, graph
