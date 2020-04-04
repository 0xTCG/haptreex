# HapTree-X

HapTree-X is a computational tool that phasees various kinds of next-generation sequencing data. 
Currently, it supports whole-genome, whole-exome, 10X Genomics and RNA-seq data. 
It is especially powerful on RNA-seq data as it can utilize allelic imbalance to better phase genic regions.


## Installation

HapTree-X binaries are available for Linux and macOS.
To download them, please go to the Releases tab and download the latest version.

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

If you want to phase 10X genomics sampes, pass `--10x` flag to HapTree-X.

Finally, HapTree-X can be run in multi-threaded mode. To enable it, set the `OMP_NUM_THREADS` variable to the desired number of threads.
An example would be:
```
OMP_NUM_THREADS=4 haptreex -v ...
```

## Paper data

Exterimental notebook and scripts used to generate the relevant paper data are located in [scripts/](scirpts) directory.

## Contact

For questions or issues, either open GitHub issue or contact us at:

Ibrahim Numanagić (inumanag@uvic.ca)
Lillian Zhang (lillianz@mit.edu)

