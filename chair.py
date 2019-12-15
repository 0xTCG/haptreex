"""
This tool extracts reads from a SAM file that overlap heterozygous-SNPs within a given VCF file and output the extracted overlap information in a fragment matrix
Fragmat output format -- Only needed for internal processing
N [read name] ([Snp No] [alleles])^N [qualities]

Default mode is going to run multi-allelic version
Biallelic mode is going to make 0/2 and print it as if it is 0 and 1 but use the second allele
"""

import sys
from dataclasses import dataclass
from typing import Tuple, Dict, List, Set, NamedTuple, Optional
import pysam
import bisect
import sys
from pprint import pprint










from collections import defaultdict
x = defaultdict(lambda: defaultdict(int))
for n, snps in process(sys.argv[1], sys.argv[2]):
    for s, a, _ in snps:
        x[s][a] += 1
for s, i in x.items():
    if len(i) < 2:
        print(f'{s} -> {dict(i)}')
    if 23504325 <= s.pos <= 23504326:
        print(f'-- {s} -> {dict(i)}')
