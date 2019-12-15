import pysam
import sys
from typing import Tuple, Dict, List, Set, NamedTuple, Optional

class SNP(NamedTuple):
    chr: str
    pos: int
    name: str
    alleles: List[str]


sample=None
vcf_path=sys.argv[1]
snps: Dict[str, List[SNP]] = {}
with pysam.VariantFile(vcf_path) as vcf:
    samples = list(vcf.header.samples)
    if not samples:
        raise ValueError('No samples present in the VCF file')
    if sample and sample not in samples:
        raise ValueError(f'Sample {sample} not found in the VCF')
    sample = sample if sample else samples[0]
    print(f'Using sample {sample} in {vcf_path}...')
    for rec in vcf:
        # We only deal with SNPs here for now
        if len(rec.ref) != 1:
            continue
        alts = [a for a in rec.alts if len(a) == 1 and a != rec.ref]

        if rec.chrom not in snps:
            snps[rec.chrom] = []

        gt = rec.samples[sample]['GT']
        # Get only alleles that are specified in GT field
        alleles = [a for i, a in enumerate(rec.alleles) if i in gt]
        if len(alleles) > 1: # Ignore homozygous SNPs
            snps[rec.chrom].append(SNP(rec.chrom, rec.pos - 1, rec.id, alleles))
for l in snps.values():
    l.sort()
return snps