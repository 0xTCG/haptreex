from gene import Gene
from read import SNP, Read
from rna import RNAData
from graph import Data, edges_from_readlist
from common import QUALITY_CUTOFF
from mytime import timing
from dataclasses import dataclass
from typing import Tuple, Dict, List, Set, NamedTuple, Optional, Iterator
import pysam
import bisect
import sys
from pprint import pprint


def parse_vcf(
    vcf_path: str,
    sample: Optional[str] = None
) -> Tuple[Dict[str, List[SNP]], Dict[int, Tuple[str, int]]]:
    """
    Parse a VCF file and return a dictionary of sorted SNPs:
        d := {chromosome_name: sorted([snp_1, snp_2])}
    and the corresponding line index
        {line_in_vcf: (chr, index in d[chr])}.
    Each snp_i is a SNP that is heterozygous in the sample (i.e. |set(GT(snp_i))| > 1).
    """

    snps: Dict[str, List[SNP]] = {}
    index: List[Tuple[str, int]] = {}
    with pysam.VariantFile(vcf_path) as vcf:
        samples = list(vcf.header.samples)
        if not samples:
            raise ValueError('No samples present in the VCF file')
        if sample and sample not in samples:
            raise ValueError(f'Sample {sample} not found in the VCF')
        sample = sample if sample else samples[0]
        print(f'Using sample {sample} in {vcf_path}...')
        for rec_i, rec in enumerate(vcf):
            if len(rec.ref) != 1: continue # Ignore indels
            # We only deal with SNPs here for now
            gt = rec.samples[sample]['GT']
            # Get only alleles that are specified in GT field
            alleles = [a for i, a in enumerate(rec.alleles) if i in gt and len(a) == 1]
            if len(alleles) > 1: # Ignore homozygous SNPs
                snp = SNP(rec_i, rec.chrom, rec.pos - 1, rec.id, alleles)
                if rec.chrom not in snps:
                    snps[rec.chrom] = [snp]
                elif snp < snps[rec.chrom][-1]:
                    raise ValueError(f'VCF is not sorted (SNP {snp})')
                else:
                    snps[rec.chrom].append(snp)
                index[rec_i] = rec.chrom, len(snps[rec.chrom])
    return snps, index


def parse_read(
    lines, #: List[SAMRecord],
    snps: Dict[str, List[SNP]],
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

    cov: Dict[SNP, Dict[int, List[str]]] = {} # SNP: {allele: [qual1, qual2, ...]}
    name = lines[0].query_name
    counts = [0] * len(lines)
    for line_i, sam in enumerate(lines):
        read, ref = 0, 0
        cand = snps[sam.reference_name]
        for op, sz in sam.cigartuples:
            if op in [0, 7, 8]: #'M=X':
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
            elif op in [1, 4]: #'IS':
                read += sz
            elif op in [2, 3, 5, 6]: #'DNHP':
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
    snps: Dict[str, List[SNP]],
    threshold: int = 1,
    ignore_chimeric: bool = True,
    ignore_duplicates: bool = True,
    ignore_conflicts: bool = True
) -> Iterator[Tuple[str, List[Tuple[SNP, int, str]]]]:
    """
    Reads a sorted SAM/BAM.
    """

    seen = {} #S Dict[str, SAMRecord]
    seen_chrs: Set[str] = {}
    with pysam.AlignmentFile(sam_path) as sam:
        for line in sam:
            if line.reference_name not in seen_chrs:
                print(f'Parsing {line.reference_name}, {len(seen)} cached so far...')
            seen_chrs.add(line.reference_name)
            name = sam.query_name
            if len(name) > 2 and name[-2] in '#/':
                name = name[:-2]

            if (
                line.is_supplementary
                or (ignore_duplicates and line.is_duplicate)
                or (ignore_chimeric and line.reference_name != line.next_reference_name)
                or line.is_unmapped
                or not line.cigartuples
                or line.reference_name not in snps
            ):
                if not line.is_supplementary and not line.mate_is_unmapped and name in seen:
                    yield from parse_read([seen[name]], snps, threshold, ignore_conflicts)
                    del seen[name]
                continue
            elif name in seen:
                yield from parse_read([seen[name], line], snps, threshold, ignore_conflicts)
                del seen[name]
            elif (
                not line.mate_is_unmapped
                and line.reference_name == line.next_reference_name
                and line.mpos < line.pos
            ):
                yield from parse_read([line], snps, threshold, ignore_conflicts)
            elif (
                not line.mate_is_unmapped
                and line.reference_name != line.next_reference_name
                and line.next_reference_name in seen_chrs
            ):
                yield from parse_read([line], snps, threshold, ignore_conflicts)
            else:
                seen[name] = line
    for line in seen:
        yield from parse_sam_pair([line], snps, threshold, ignore_conflicts)


def parse_fragmat(
    fragmat: str,
    snp_index: Dict[int, Tuple[str, int]],
    snps: Dict[str, List[SNP]],
    skip_single: bool
) -> Iterator[Tuple[str, List[Tuple[SNP, int, str]]]]:
    # translates fragment matrix into a List of reads
    print(f"Loading and formatting fragments...")
    read_list_list: List[List[Tuple[int, int]]] = []
    with open(fragmat, "r") as f:
        for r in f:
            frag = r.split()
            name, qual, frag = frag[1] frag[-1], frag[2:-1]
            if len(frag) % 2 == 1:
                raise ValueError(f"fragment file error: {frag}")
            if len(qual) == 1 and ord(qual) < QUALITY_CUTOFF: # TODO: fix this
                continue

            alleles = []
            for i in range(0, len(frag), 2):
                if int(frag[i]) not in snp_index:
                    raise ValueError(f'Invalid SNP index {frag[i]} for read {name}')
                chr, idx = snp_index[int(frag[i])]
                for j, allele in enumerate(frag[i + 1]):
                    s = snps[chr][idx + j - 1]
                    if not 0 <= int(allele) < len(s.alleles):
                        raise ValueError(f'Invalid allele {allele} for SNP {s}')
                    alleles.append((s, int(allele), qual[len(alleles)]))
            yield name, alleles


def parse_phases(
    vcf: Tuple[Dict[str, List[SNP]], Dict[int, Tuple[str, int]]],
    paths: List[str],
    skip_single: bool = True
) -> Iterator[Read]:
    snps, snp_index = vcf
    reads = []
    for path in paths:
        print(f'Parsing {path}...')
        for _, alleles in parse_fragmat(path, snp_index, snps, skip_single):
            if not (skip_single and len(alleles) <= 1):
                reads.append(alleles)
    print(f"{len(reads)} reads of sufficient quality")

    read_counter = {} #S : Dict[List[Tuple[int, int]], int] = {}
    for r in reads: #S
        key = tuple((s, a) for s, a, _ in r)
        if key not in read_counter:
            read_counter[key] = 1
        else:
            read_counter[key] += 1
    print(f"{len(read_counter)} distinct reads")

    for i, (tup, cnt) in enumerate(read_counter.items()):
        yield Read(dict(tup), cnt, i)


###########################################################################
# make RNA_data and DNA_data objects


# def make_RNA_data_from_fragmat(
#     gene_data: str,
#     fragmats: List[str],
#     vcf: str,
#     error: float,
#     isoforms: str
# ) -> RNAData:
#     read_list = make_readlist_from_fragmat(fragmats, skip_single=False) # 4s
#     print(f"Loading VCF file {vcf}")
#     S, names, chroms, positions, k = positions_names_states(vcf) # 3s
#     print("Preparing data for ReadGraph")
#     genes = determine_genes_gtf(gene_data, set(chroms.values())) # 42s

#     isodict: Dict[str, Tuple[float, str]] = {}
#     filtered_genes = genes
#     if isoforms != "":
#         print("Building IsoDict")
#         isodict = build_isodict(isoforms)
#         filtered_genes = filter_transcripts(genes, isodict)

#     return RNAData(
#         S, genes, filtered_genes, error, read_list, positions, names, chroms, isodict
#     )


def load_dna_data(
    vcf_path: str,
    paths: List[str],
    error: float,
    RNA_readlist: Dict[int, Read] = None
) -> Graph:
    print(f"Loading VCF file {vcf_path}"...)
    vcf = parse_vcf(vcf_path)
    print(f"{len(vcf[1])} SNPs in VCF file")

    reads = list(parse_phases(vcf, paths, skip_single=True))
    # if len(RNA_readlist) > 0:
    #     max_key = max(reads.keys())
    #     for zz in RNA_readlist:
    #         reads[zz + max_key] = RNA_readlist[zz]
    for r in reads.values():
        r.special_snp = sorted(r.snps)[1]
        r.rates = [0.5, 0.5]

    return Graph(reads, ploidy=2, error=error)
