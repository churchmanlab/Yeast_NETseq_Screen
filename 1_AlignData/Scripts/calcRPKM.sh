#!/bin/bash

# Script to calculate gene-by-gene RPKM for given deletion strain

# Written by: K. Lachance
# Date: November 26, 2018

# Use: ./calcRPKM.sh mut

# Parameters
mut=$1

# Make output directory
mkdir -p RPKMnormFiles

# Use bedtools to calculate RPKM for each strand
cat ./RPMnormFiles/${mut}.pos.bed | grep -v track | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' | sortBed -i - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding_pos.bed - | awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7/(($3-$2)/1000)}' > ./RPKMnormFiles/${mut}.RPKMnorm.pos.bed

cat ./RPMnormFiles/${mut}.neg.bed | grep -v track | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' | sortBed -i - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding_neg.bed - | awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7/(($3-$2)/1000)}' > ./RPKMnormFiles/${mut}.RPKMnorm.neg.bed

# Combine strands and sort
cat ./RPKMnormFiles/${mut}.RPKMnorm.pos.bed ./RPKMnormFiles/${mut}.RPKMnorm.neg.bed | sortBed -i - > ./RPKMnormFiles/${mut}.RPKMnorm.bed
