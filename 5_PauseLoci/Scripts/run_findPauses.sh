#!/bin/bash

# Script to identify pause locations in genes using NET-seq data (as defined by Churchman & Weissman, 2011)

# Written by: K. Lachance
# Date: November 1, 2018

# Use: ./run_findPauses.sh mut

# Mutant of interest
mut=$1

# Echo deletion strain
echo $mut

# Get bedgraph files
posSignal=/n/groups/churchman/kcl19/Screen/DATA/NonormBedgraphFiles/${mut}-1.pos.bed # For replicates
negSignal=/n/groups/churchman/kcl19/Screen/DATA/NonormBedgraphFiles/${mut}-1.neg.bed

# Get gene files
posGenes=../0_Annotations/Genes/genes_all_pos.bed
negGenes=../0_Annotations/Genes/genes_all_neg.bed

# Overlap and isolate genes that are above coverage threshold
t=5
awk -F "\t" 'BEGIN {OFS="\t"} ($4 > 0) {print $1, $2, $3, $4*($3-$2)}' $posSignal | mapBed -a $posGenes -b - -c 4 -o sum | awk -F "\t" 'BEGIN {OFS="\t"} ($7 != ".") {print $0, $7/($3-$2)}' | awk -v t=$t 'BEGIN {OFS="\t"} ($8 > t) {print $1, $2, $3, $4, $5, $6}' > ${mut}-1.covThreshold${t}.pos.bed
awk -F "\t" 'BEGIN {OFS="\t"} ($4 > 0) {print $1, $2, $3, $4*($3-$2)}' $negSignal | mapBed -a $negGenes -b - -c 4 -o sum | awk -F "\t" 'BEGIN {OFS="\t"} ($7 != ".") {print $0, $7/($3-$2)}' | awk -v t=$t 'BEGIN {OFS="\t"} ($8 > t) {print $1, $2, $3, $4, $5, $6}' > ${mut}-1.covThreshold${t}.neg.bed

# Repeat for replicate 2
posSignal=/n/groups/churchman/kcl19/Screen/DATA/NonormBedgraphFiles/${mut}-2.pos.bed
negSignal=/n/groups/churchman/kcl19/Screen/DATA/NonormBedgraphFiles/${mut}-2.neg.bed

awk -F "\t" 'BEGIN {OFS="\t"} ($4 > 0) {print $1, $2, $3, $4*($3-$2)}' $posSignal | mapBed -a $posGenes -b - -c 4 -o sum | awk -F "\t" 'BEGIN {OFS="\t"} ($7 != ".") {print $0, $7/($3-$2)}' | awk -v t=$t 'BEGIN {OFS="\t"} ($8 > t) {print $1, $2, $3, $4, $5, $6}' > ${mut}-2.covThreshold${t}.pos.bed
awk -F "\t" 'BEGIN {OFS="\t"} ($4 > 0) {print $1, $2, $3, $4*($3-$2)}' $negSignal | mapBed -a $negGenes -b - -c 4 -o sum | awk -F "\t" 'BEGIN {OFS="\t"} ($7 != ".") {print $0, $7/($3-$2)}' | awk -v t=$t 'BEGIN {OFS="\t"} ($8 > t) {print $1, $2, $3, $4, $5, $6}' > ${mut}-2.covThreshold${t}.neg.bed

# Only take overlapping genes from two replicates
grep -f ${mut}-1.covThreshold${t}.pos.bed ${mut}-2.covThreshold${t}.pos.bed > ${mut}.covThreshold${t}.pos.bed
grep -f ${mut}-1.covThreshold${t}.neg.bed ${mut}-2.covThreshold${t}.neg.bed > ${mut}.covThreshold${t}.neg.bed

# Redefine gene files
posGenes=${mut}.covThreshold${t}.pos.bed
negGenes=${mut}.covThreshold${t}.neg.bed

# Counter for progress
Npos=`wc -l $posGenes | awk '{print $1}'`
i=1

# For each gene (pos)
while read gene
do
    # Progress tracker
    if (( $i % 100 == 0 )) 
    then
	echo "${i} / ${Npos} genes done!"
    fi

    # Store gene as one-line to temporary file
    echo $gene > ${mut}.tmpGene.bed

    # Overlap gene with NET-seq signal, with one row per basepair
    bedops --chop ${mut}.tmpGene.bed | intersectBed -a $posSignal -b - > ${mut}.tmpCov.bed

    # Call R script to find pauses
    ./Scripts/findPauses_NB.R ${mut}.tmpCov.bed ${mut} 3 pos

    #  Increment tracker
    i=$((i+1))

done < $posGenes
echo "Done with genes on the + strand!"

# Same for genes on neg strand
Nneg=`wc -l $negGenes | awk '{print $1}'`
i=1

while read gene
do
    if (( $i % 100 == 0 ))
    then
        echo "${i} / ${Nneg} genes done!"
    fi

    echo $gene > ${mut}.tmpGene.bed
    bedops --chop ${mut}.tmpGene.bed | intersectBed -a $negSignal -b - > ${mut}.tmpCov.bed

    ./Scripts/findPauses_NB.R ${mut}.tmpCov.bed ${mut} 3 neg

    i=$((i+1))

done < $negGenes
echo "Done with genes on the - strand!"

# Remove temporary files
#rm -fr tmp???.bed
