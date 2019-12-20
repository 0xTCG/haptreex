from read import SNP, Read
from graph import Graph
from dataclasses import dataclass
from rna import Gene, RNAGraph
import os
import pysam
from typing import Tuple, Dict, List, Set, Optional, Iterator, Any


QUALITY_CUTOFF = 10


@dataclass
class VCF:
    snps: List[SNP]
    line_to_snp: Dict[int, int]  # 1-based index
    chromosomes: Dict[str, Tuple[int, int]]  # Chromosome name to the SNP span

    def find_first(self, chr: str, pos: int) -> int:
        if chr not in self.chromosomes:
            return -1
        lo, hi = self.chromosomes[chr]
        while lo < hi:
            mid = (lo + hi) // 2
            if self.snps[mid].pos < pos:
                lo = mid + 1
            else:
                hi = mid
        return lo


def parse_vcf(vcf_path: str, sample: Optional[str] = None) -> VCF:
    """
    Parse a VCF file and return a dictionary of sorted SNPs:
        d := {chromosome_name: sorted([snp_1, snp_2])}
    and the corresponding line index
        {line_in_vcf: (chr, index in d[chr])}.
    Each snp_i is a SNP that is heterozygous in the sample (i.e. |set(GT(snp_i))| > 1).
    """

    snps: List[SNP] = []
    line_to_snp: Dict[int, int] = {}
    chromosomes: Dict[str, Tuple[int, int]] = {}
    with open(vcf_path) as vcf:
        # samples = list(vcf.header.samples)
        # if not samples:
        #     raise ValueError('No samples present in the VCF file')
        # if sample and sample not in samples:
        #     raise ValueError(f'Sample {sample} not found in the VCF')
        # sample = sample if sample else samples[0]
        print(f'Using sample {sample} in {vcf_path}...')
        prev_chr, seen_snps = "", 0
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
                line_to_snp[seen_snps] = len(snps)
                if chr not in chromosomes:
                    if prev_chr:
                        chromosomes[prev_chr] = (chromosomes[prev_chr][0], len(snps))
                    chromosomes[chr] = (len(snps), -1)
                    prev_chr = chr
                snps.append(snp)
    if prev_chr:
        chromosomes[prev_chr] = (chromosomes[prev_chr][0], len(snps))
    return VCF(snps, line_to_snp, chromosomes)


def parse_gtf(gtf_path: str, chroms: Dict[str, Tuple[int, int]]) -> Iterator[Gene]:
    b = [a[:-2].split("\t") for a in open(gtf_path, "r") if a[0] != '#']

    def parse_f(f):
        x, y = f.split(" ", 1)
        return x, y[1:-1] if y[0] == y[-1] == '"' else y

    gi, i = 0, 0
    while i < len(b):
        if b[i][2] == "transcript":
            chr, interval, sign, _info = (
                b[i][0], (int(b[i][3]) - 1, int(b[i][4])), b[i][6], b[i][8]
            )
            if chr not in chroms:
                continue
            info = dict(parse_f(f) for f in _info.split("; "))
            name = info.get("gene_name", info["transcript_id"])
            i += 1
            exons: List[Tuple[int, int]] = []
            while i < len(b) and b[i][2] != "transcript":
                if b[i][2] == "exon":
                    exon_start, exon_end = int(b[i][3]) - 1, int(b[i][4])
                    exons.append((exon_start, exon_end))
                i += 1
            # print('Adding', gi, chr, info['gene_id'], name, len(exons))
            yield Gene(gi, name, chr, interval, sign, exons)
            gi += 1
        else:
            i += 1


def parse_read(
    vcf: VCF,
    lines,
    threshold: int,
    ignore_conflicts: bool = True
) -> Iterator[Tuple[str, List[Tuple[int, int, str]]]]:
    """
    If reads are valid and pass the threshold filter, yields the
        (read_name, [all_1, all_2, ...])
    where
        all_i := (snp, allele, quality).
    Example:
        ('read1', [(SNP("chr1", 12), 'A', 'E'), ...])
    """

    cov: Dict[int, Dict[int, List[str]]] = {}  # SNP: {allele: [qual1, qual2, ...]}
    name = lines[0].query_name
    counts = [0] * len(lines)
    for line_i, sam in enumerate(lines):
        read, ref = 0, 0
        for op, sz in sam.cigartuples:
            if op in [0, 7, 8]:  # 'M=X':
                start = vcf.find_first(sam.reference_name, sam.pos + ref)
                for i in range(start, vcf.chromosomes[sam.reference_name][1]):
                    snp = vcf.snps[i]
                    if snp.pos >= sam.pos + ref + sz:
                        break
                    t = snp.pos - sam.pos - ref + read
                    if sam.seq[t] in snp.alleles:
                        allele = snp.alleles.index(sam.seq[t])
                        qual = sam.qual[t] if sam.qual != '*' else '.'
                        cov.setdefault(i, {}).setdefault(allele, []).append(qual)
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
    if ignore_conflicts:
        for i in list(cov):
            if len(cov[i]) > 1:
                del cov[i]
    if len(cov) >= threshold:
        yield name, [
            (sid, a, max(q))
            for sid, A in cov.items()
            for a, q in A.items()
        ]


def parse_bam(
    vcf: VCF,
    sam_path: str,
    threshold: int = 1,
    no_duplicates: bool = True,
    no_conflicts: bool = True
) -> Iterator[Tuple[str, List[Tuple[int, int, str]]]]:
    """
    Reads a sorted SAM/BAM.
    """

    seen: Dict[str, Any] = {}  # Dict[str, SAMRecord]
    seen_chrs: Set[str] = set()
    with pysam.AlignmentFile(sam_path) as sam:
        for line_i, line in enumerate(sam):
            if line.reference_name not in seen_chrs or line_i % 100000 == 0:
                pass
                # print(
                #     f'\rParsing {line.reference_name}; '
                #     f'{len(seen):,} cached so far; '
                #     f'{line_i:,} processed so far...', end='')
            seen_chrs.add(line.reference_name)
            name = line.query_name
            if len(name) > 2 and name[-2] in '#/':
                name = name[:-2]

            if (
                line.is_supplementary
                or (no_duplicates and line.is_duplicate)
                or line.reference_name != line.next_reference_name
                or line.is_unmapped
                or not line.cigartuples
                or line.reference_name not in vcf.chromosomes
            ):
                if not line.is_supplementary and not line.mate_is_unmapped and name in seen:
                    yield from parse_read(vcf, [seen[name]], threshold, no_conflicts)
                    del seen[name]
                continue
            elif name in seen:
                yield from parse_read(vcf, [seen[name], line], threshold, no_conflicts)
                del seen[name]
            # elif (
            #     not line.mate_is_unmapped
            #     and line.reference_name == line.next_reference_name
            #     and line.mpos < line.pos
            # ):
            #     yield from parse_read(vcf, [line], threshold, no_conflicts)
            else:
                seen[name] = line
    for line in seen.values():
        yield from parse_read(vcf, [line], threshold, no_conflicts)


def parse_fragmat(
    vcf: VCF, fragmat: str
) -> Iterator[Tuple[str, List[Tuple[int, int, str]]]]:
    """
    Yields reads in the format
        (read_name, [(snp_id, allele_id, qual), ...])
    """
    print(f"Loading and formatting fragments...")
    with open(fragmat, "rb") as f:
        header = f.read(1024)
        textchars = {7, 8, 9, 10, 12, 13, 27} | set(range(0x20, 0x100)) - {0x7f}
        if any(c not in textchars for c in header):
            raise ValueError(f'{fragmat} is a binary file')
        if not chr(header[0]).isdigit():
            raise ValueError(f'{fragmat} is not valid fragmat file')
    with open(fragmat) as f:
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
                    alleles.append((snp.id, int(allele), qual[len(alleles)]))
                    idx += 1
            yield name, alleles


def parse_phases(
    vcf: VCF,
    paths: List[str],
    skip_single: bool = True
) -> Iterator[Read]:
    reads = []
    coverages = {}
    for path in paths:
        print(f'Parsing {path}...')
        _, ext = os.path.splitext(path)
        if ext in ['.sam', '.bam', '.cram']:
            it = parse_bam(vcf, path)
        else:
            it = parse_fragmat(vcf, path)
        for name, alleles in it:
            if not (skip_single and len(alleles) <= 1):
                for pos, al, _ in alleles:
                    coverages.setdefault(pos, set()).add(al)
                reads.append(alleles)
    print(f"{len(reads)} reads of sufficient quality")

    read_counter = {}  #S : Dict[List[Tuple[int, int]], int] = {}
    for r in reads:
        key = tuple(sorted(
            (s, a) for s, a, _ in r if s in coverages and len(coverages[s]) >= 2
        ))
        if len(key) <= int(skip_single):
            continue
        if key not in read_counter:
            read_counter[key] = 1
        else:
            read_counter[key] += 1
    print(f"{len(read_counter)} distinct reads")

    for i, (tup, cnt) in enumerate(read_counter.items()):
        yield Read({snp_id: snp_a for snp_id, snp_a in tup}, cnt, i)


def load_rna_data(
    vcf: VCF,
    gtf_path: str,
    paths: List[str],
    isoforms_path: str
) -> RNAGraph:
    print(f"Loading GTF {gtf_path}...")
    genes = list(parse_gtf(gtf_path, vcf.chromosomes))
    print(f"{len(genes)} genes in GTF file")

    reads = list(parse_phases(vcf, paths, skip_single=False))
    if isoforms_path:
        print("Building IsoDict...")
        isodict = build_isodict(isoforms_path)
        genes = filter_transcripts(genes, isodict)

    return RNAGraph(vcf, genes, reads)


def load_dna_data(
    vcf: VCF,
    paths: List[str],
    rna_reads: List[Read] = None
) -> Graph:
    reads = list(parse_phases(vcf, paths, skip_single=True))
    if rna_reads:
        for rna_read in rna_reads:
            reads.append(rna_read)
    for r in reads:
        r.special_snp = sorted(r.snps)[1]
        r.rates = [0.5, 0.5]

    return Graph(reads, ploidy=2)
