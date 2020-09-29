#!/bin/bash

# Script calculates the AS to S ratio for gene-body and tandem regions of the genome

# Written by: K. Lachance
# Date: October 29, 2018

# Use: ./calcAStoSratio.sh mut

# Input deletion strain name
mut=$1

# Echo to terminal
echo "Working on ${mut} deletion strain..."

# Signal files with NET-seq data
# TODO: When newly aligned bedfiles are created, point to those instead of the old ones (now)
posSigFile=/n/groups/churchman/kcl19/Screen/DATA/CombinedReplicates/${mut}.pos.bed
negSigFile=/n/groups/churchman/kcl19/Screen/DATA/CombinedReplicates/${mut}.neg.bed

# Annotation files for genes (coding and tandem)
# Note: "TandemUp" is defined as the region -100 to -500 upstream of the TSS and "TandemDown" is the region 0 to +500 from the TSS
posGeneFile=../0_Annotations/Genes/genes_polII.noOver.coding_pos.bed
negGeneFile=../0_Annotations/Genes/genes_polII.noOver.coding_neg.bed

posTandemUpFile=../0_Annotations/Genes/genes_polII.noOver.coding.TandemUp_pos.bed
negTandemUpFile=../0_Annotations/Genes/genes_polII.noOver.coding.TandemUp_neg.bed
posTandemDownFile=../0_Annotations/Genes/genes_polII.noOver.coding.TandemDown_pos.bed
negTandemDownFile=../0_Annotations/Genes/genes_polII.noOver.coding.TandemDown_neg.bed

# Make directories if needed
mkdir -p CodingAntisense
mkdir -p TandemAntisense

# Overlap signal and annotation files
# Gene, + strand
cat $posSigFile | awk -F "\t" '{print $1, $2, $3, ".", $4}' | sed 's/ /\t/g' | bedmap --sum $posGeneFile - > tmp.gene_sense.txt
cat $negSigFile | awk -F "\t" '{print $1, $2, $3, ".", $4}' | sed 's/ /\t/g' | bedmap --sum $posGeneFile - > tmp.gene_anti.txt
paste $posGeneFile tmp.gene_sense.txt tmp.gene_anti.txt | awk -F "\t" '{if ($7 + 0 != 0) print $0, $8/$7; else print $0, "NA"}' |sed 's/ /\t/g' > CodingAntisense/$mut.geneNETseqRatio.pos.bed # File has 3 additional columns: sense, antisense, and AS/S ratio

# Tandem, + strand
cat $posSigFile | awk -F "\t" '{print $1, $2, $3, ".", $4}' | sed 's/ /\t/g' | bedmap --sum $posTandemDownFile - > tmp.tandem_sense.txt
cat $negSigFile | awk -F "\t" '{print $1, $2, $3, ".", $4}' | sed 's/ /\t/g' | bedmap --sum $posTandemUpFile - > tmp.tandem_anti.txt
paste $posTandemDownFile tmp.tandem_sense.txt tmp.tandem_anti.txt | awk -F "\t" '{if ($7 + 0 != 0) print $0, $8/$7; else print $0, "NA"}' |sed 's/ /\t/g' > TandemAntisense/$mut.tandemNETseqRatio.pos.bed


# Genes, - strand
cat $negSigFile | awk -F "\t" '{print $1, $2, $3, ".", $4}' | sed 's/ /\t/g' | bedmap --sum $negGeneFile - > tmp.gene_sense.txt
cat $posSigFile | awk -F "\t" '{print $1, $2, $3, ".", $4}' | sed 's/ /\t/g' | bedmap --sum $negGeneFile - > tmp.gene_anti.txt
paste $negGeneFile tmp.gene_sense.txt tmp.gene_anti.txt | awk -F "\t" '{if ($7 + 0 != 0) print $0, $8/$7; else print $0, "NA"}' |sed 's/ /\t/g' > CodingAntisense/$mut.geneNETseqRatio.neg.bed

# Tandem, - strand
cat $negSigFile | awk -F "\t" '{print $1, $2, $3, ".", $4}' | sed 's/ /\t/g' | bedmap --sum $negTandemDownFile - > tmp.tandem_sense.txt
cat $posSigFile | awk -F "\t" '{print $1, $2, $3, ".", $4}' | sed 's/ /\t/g' | bedmap --sum $negTandemUpFile - > tmp.tandem_anti.txt
paste $negTandemDownFile tmp.tandem_sense.txt tmp.tandem_anti.txt | awk -F "\t" '{if ($7 + 0 != 0) print $0, $8/$7; else print $0, "NA"}' |sed 's/ /\t/g' > TandemAntisense/$mut.tandemNETseqRatio.neg.bed

# Cleanup
rm -fr tmp.*.txt
