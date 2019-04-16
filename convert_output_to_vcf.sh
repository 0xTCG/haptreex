#!/bin/bash

# ARG 1 -> HapTree-X output
# ARG 2 -> VCF used as input to HapTree-X
# ARG 3 -> Output VCF that contains Haptree-X blocks 

awk 'BEGIN{block_index = 0} NR == FNR { if($1 ~ "BLOCK") {block_index++}; if($1 !~ "BLOCK" && NF == 5) {block[$1]=block_index; phase[$1]=$2; chrom[$1]=$4; pos[$1]=$5;}} NR != FNR { if(substr($1,1,2) == "##") {print $0} if(substr($1,1,1) == "#" && substr($1,2,1) != "#") {print "##FORMAT=<ID=HXB,Number=1,Type=Integer,Description=\"HapTree-X Phased Haplotype Block Index\">"; print "##FORMAT=<ID=HXP,Number=1,Type=String,Description=\"HapTree-X Phased Genotype Within Block\">"; print $0; hx_index = 0;} if(substr($1,1,1) != "#") {hx_index++; if(length(chrom[hx_index]) > 0) {if(chrom[hx_index] == $1 && pos[hx_index] == $2) {HXB=block[hx_index]; HXP=phase[hx_index] "|" 1-phase[hx_index];} else {print "ERROR: HapTree-X output is inconsistent with VCF line at index " hx_index " (printed below), make sure the VCF input file used for HapTree-X and this script are the same.\n" $0 > "/dev/stderr"; exit 1;}} for(i=1; i<NF-1; i++) {printf $i "\t"} printf $(NF-1); if(HXB!=0) {printf ":HXB:HXP"} printf "\t"; printf $NF; if(HXB!=0) {printf ":" HXB ":" HXP} printf "\n"; HXB = 0; HXP = ""}}' $1 $2 > $3
