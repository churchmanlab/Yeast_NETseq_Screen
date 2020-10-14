#!/bin/bash

# Calculate % of reads in highly-expressed genes in pauses
# Use: ./calcPauseStrength.sh deletionStrains.txt

while read mut
do

    echo ${mut}

    # Get bedgraph files
    posSignal=../1_AlignData/CoverageFiles/${mut}-1_stt_pos.bedGraph 
    negSignal=../1_AlignData/CoverageFiles/${mut}-1_stt_neg.bedGraph 

    # Get gene files
    posGenes=${mut}.covThreshold5.pos.bed
    negGenes=${mut}.covThreshold5.neg.bed

    # Overlap and isolate genes that are above coverage threshold
    awk -F "\t" 'BEGIN {OFS="\t"} ($4 > 0) {print $1, $2, $3, $4*($3-$2)}' $posSignal | mapBed -a $posGenes -b - -c 4 -o sum | awk -F "\t" 'BEGIN {OFS="\t"} ($7 != ".") {print $0, $7/($3-$2)}' | awk -v t=$t 'BEGIN {OFS="\t"} ($8 > t) {print $1, $2, $3, $4, $5, $6, $7, $8}' > ${mut}-1.covThreshold${t}.pos.bed
    awk -F "\t" 'BEGIN {OFS="\t"} ($4 > 0) {print $1, $2, $3, $4*($3-$2)}' $negSignal | mapBed -a $negGenes -b - -c 4 -o sum | awk -F "\t" 'BEGIN {OFS="\t"} ($7 != ".") {print $0, $7/($3-$2)}' | awk -v t=$t 'BEGIN {OFS="\t"} ($8 > t) {print $1, $2, $3, $4, $5, $6, $7, $8}' > ${mut}-1.covThreshold${t}.neg.bed

    # Repeat for replicate 2
    posSignal=../1_AlignData/CoverageFiles/${mut}-2_stt_pos.bedGraph 
    negSignal=../1_AlignData/CoverageFiles/${mut}-2_stt_neg.bedGraph 

    awk -F "\t" 'BEGIN {OFS="\t"} ($4 > 0) {print $1, $2, $3, $4*($3-$2)}' $posSignal | mapBed -a $posGenes -b - -c 4 -o sum | awk -F "\t" 'BEGIN {OFS="\t"} ($7 != ".") {print $0, $7/($3-$2)}' | awk -v t=$t 'BEGIN {OFS="\t"} ($8 > t) {print $1, $2, $3, $4, $5, $6, $7, $8}' > ${mut}-2.covThreshold${t}.pos.bed
    awk -F "\t" 'BEGIN {OFS="\t"} ($4 > 0) {print $1, $2, $3, $4*($3-$2)}' $negSignal | mapBed -a $negGenes -b - -c 4 -o sum | awk -F "\t" 'BEGIN {OFS="\t"} ($7 != ".") {print $0, $7/($3-$2)}' | awk -v t=$t 'BEGIN {OFS="\t"} ($8 > t) {print $1, $2, $3, $4, $5, $6, $7, $8}' > ${mut}-2.covThreshold${t}.neg.bed

    # Count number of reads in possible pauses for both replicates
    reads1=`cat ${mut}-1.covThreshold5.???.bed | awk -F "\t" '{sum += $7} END {print sum}'`
    reads2=`cat ${mut}-3.covThreshold5.???.bed | awk -F "\t" '{sum += $7} END {print sum}'`

    # Count number of reads in pauses for both replicates
    pauses1=`cat ../../${mut}_PASS_pauseScore.txt | awk -F "\t" '{sum += $4} END {print sum}'`
    pauses2=`cat ../../${mut}_PASS_pauseScore.txt | awk -F "\t" '{sum += $5} END {print sum}'`

    echo -e "${mut}\t${reads1}\t${reads2}\t${pauses1}\t${pauses2}" >> pauseStrength.txt

done < $1

./Scripts/plotPauseStrength.R
