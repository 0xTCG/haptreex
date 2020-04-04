# HapTree-X

HapTree-X is a computational tool that phases various kinds of next-generation sequencing data. 
Currently, it supports whole-genome, whole-exome, 10X Genomics and RNA-seq data. 
It is especially powerful on RNA-seq data as it can utilize allelic imbalance to better phase genic regions.

## Installation

HapTree-X binaries are available for Linux and macOS in [build](build) directory.

You will need libomp to run HapTree-X. 
`brew install libomp` (macOS) or `apt install libomp5` (Debian/Ubuntu) should do the trick.

## Building

To build HapTree-X from scratch, you will need the latest version of [Seq compiler](seq-lang.org) and [llc](https://llvm.org/docs/CommandGuide/llc.html) (typically shipped with LLVM).

Then, issue `make` to build HapTree-X. The resulting binary will be located in `build/haptreex`.

## Usage

Basic usage is:
```
haptreex -v [VCF file with variants to be phased]
         -d [indexed SAM, BAM or CRAM file with the aligned reads]
         -o [output file]
```

For RNA-seq data, use:
```
haptreex -v [VCF file with variants to be phased]
         -r [indexed SAM, BAM or CRAM file with the aligned reads]
         -g [GTF file compatible with the provided SAM/BAM/CRAM]
         -o [output file]
```
`-g` parameter is optional. However, its inclusion will result in better phases.

HapTree-X can also phase both DNA and RNA-seq samples at the same time. An example would be:
```
haptreex -v [VCF file with variants to be phased]
         -r [indexed RNA-seq SAM, BAM or CRAM file with the aligned reads]
         -d [indexed WGS/WXS SAM, BAM or CRAM file with the aligned reads]
         -g [GTF file compatible with the provided SAM/BAM/CRAM]
         -o [output file]
```

If you want to phase 10X genomics samples, pass `--10x` flag to HapTree-X.

Finally, HapTree-X can be run in multi-threaded mode. To enable it, set the `OMP_NUM_THREADS` variable to the desired number of threads.
An example would be:
```
OMP_NUM_THREADS=4 haptreex -v ...
```

## Output

HapTree-X's output follows the [HapCUT](https://github.com/vibansal/HapCUT2) output format convention. The output file will contain the set of phased haplotype blocks in a list format where the beginning of each block starts with `BLOCK` and the end of each block is indicated by `*****`.

Each line in between contains 5 tab-delimited fields, which are in order:
1. Line number in the VCF file (ignoring header lines) that contain the het-SNP
2. Phase of the het-SNP corresponding to the first digit in 0|1 or 1|0
3. Phase of the het-SNP corresponding to the second digit in 0|1 or 1|0
4. Chromosome name
5. Chromosome position 

## Paper data

Experimental notebook and the scripts used to generate the relevant paper data are located in [scripts/](scripts) directory.

## Contact

For questions or issues, either open GitHub issue or contact us at:

- Ibrahim Numanagić (inumanag at uvic dot canada)
- Lillian Zhang (lillianz at mit dot education)

