#!/bin/bash

mut=$1

# For real pauses
cat ${mut}_PASS_pauseScore.txt | awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, $6}' > ${mut}.IDRrep.PAUSES.bed 
cat ${mut}_PASS_pauseScore.txt | awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2-10, $3+10, $6}' > ${mut}.IDRrep.pauseWindows.bed
bedtools getfasta -fi ../0_Annotations/Scerevisiae_R64/sacCer3.fa -bed ${mut}.IDRrep.pauseWindows.bed -s > ${mut}.IDRrep.pauseWindows.fasta

# To shuffle pauses
cat ${mut}.covThreshold5.pos.bed ${mut}.covThreshold5.neg.bed | sortBed -i > ${mut}.covThreshold5.bed
cat ${mut}.covThreshold5.bed | grep "+" > ${mut}.covThreshold5.pos.bed
cat ${mut}.covThreshold5.bed | grep -v "+"> ${mut}.covThreshold5.neg.bed
cat ${mut}.IDRrep.PAUSES.bed | grep "+" | sortBed -i - > ${mut}.IDRrep.PAUSES.pos.bed 
cat ${mut}.IDRrep.PAUSES.bed | grep -v "+" | sortBed -i - > ${mut}.IDRrep.PAUSES.neg.bed

mergeBed -i ${mut}.covThreshold5.pos.bed | intersectBed -a - -b ${mut}.IDRrep.PAUSES.pos.bed -c > ${mut}.pausesByGene.pos.bed
mergeBed -i ${mut}.covThreshold5.neg.bed | intersectBed -a - -b ${mut}.IDRrep.PAUSES.neg.bed -c > ${mut}.pausesByGene.neg.bed

./Scripts/shufflePauses.R ${mut} pos
./Scripts/shufflePauses.R ${mut} neg

cat ${mut}.IDRrep.SHUFFLE.pos.bed ${mut}.IDRrep.SHUFFLE.neg.bed | sortBed -i > ${mut}.IDRrep.SHUFFLE.bed

# Organize files
#mkdir -p ${mut}
#mv ${mut}*bed ${mut}*fasta ${mut}

