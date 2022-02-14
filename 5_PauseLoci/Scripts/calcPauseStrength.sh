#!/bin/bash

#SBATCH -c 1
#SBATCH -t 0-01:00
#SBATCH -p short
#SBATCH --mem=10G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<youremail@address.com>
###

# Calculate % of reads in highly-expressed genes in pauses
# Use: 
# t=2
# sbatch -o logs/5_calcPauseStrength.log -e logs/5_calcPauseStrength.err ./Scripts/calcPauseStrength.sh deletionStrains.txt $t

t=$2

while read mut
do

    echo ${mut}

    # Get bedgraph files
    posSignal=../1_AlignData/${mut}-1_stt_pos_sort.bedGraph 
    negSignal=../1_AlignData/${mut}-1_stt_neg_sort.bedGraph 

    # Get gene files
    posGenes=${mut}.covThreshold${t}.pos.bed
    negGenes=${mut}.covThreshold${t}.neg.bed

	# Genome file for chr ordering
	genome=../0_Annotations/Scerevisiae_R64/SacCer3_R64_genomic.fasta.fai

    # Overlap and isolate genes that are above coverage threshold
    awk -F "\t" 'BEGIN {OFS="\t"} ($4 > 0) {print $1, $2, $3, $4*($3-$2)}' $posSignal | mapBed -a $posGenes -g $genome -b - -c 4 -o sum | awk -F "\t" 'BEGIN {OFS="\t"} ($7 != ".") {print $0, $7/($3-$2)}' | awk -v t=$t -F "\t" 'BEGIN {OFS="\t"} ($8 > t) {print $1, $2, $3, $4, $5, $6, $7, $8}' > ${mut}-1.covThreshold${t}.pos.bed
    awk -F "\t" 'BEGIN {OFS="\t"} ($4 > 0) {print $1, $2, $3, $4*($3-$2)}' $negSignal | mapBed -a $negGenes -g $genome -b - -c 4 -o sum | awk -F "\t" 'BEGIN {OFS="\t"} ($7 != ".") {print $0, $7/($3-$2)}' | awk -v t=$t -F "\t" 'BEGIN {OFS="\t"} ($8 > t) {print $1, $2, $3, $4, $5, $6, $7, $8}' > ${mut}-1.covThreshold${t}.neg.bed

    # Repeat for replicate 2
    posSignal=../1_AlignData/${mut}-2_stt_pos_sort.bedGraph 
    negSignal=../1_AlignData/${mut}-2_stt_neg_sort.bedGraph 

    awk -F "\t" 'BEGIN {OFS="\t"} ($4 > 0) {print $1, $2, $3, $4*($3-$2)}' $posSignal | mapBed -a $posGenes -g $genome -b - -c 4 -o sum | awk -F "\t" 'BEGIN {OFS="\t"} ($7 != ".") {print $0, $7/($3-$2)}' | awk -v t=$t -F "\t" 'BEGIN {OFS="\t"} ($8 > t) {print $1, $2, $3, $4, $5, $6, $7, $8}' > ${mut}-2.covThreshold${t}.pos.bed
    awk -F "\t" 'BEGIN {OFS="\t"} ($4 > 0) {print $1, $2, $3, $4*($3-$2)}' $negSignal | mapBed -a $negGenes -g $genome -b - -c 4 -o sum | awk -F "\t" 'BEGIN {OFS="\t"} ($7 != ".") {print $0, $7/($3-$2)}' | awk -v t=$t -F "\t" 'BEGIN {OFS="\t"} ($8 > t) {print $1, $2, $3, $4, $5, $6, $7, $8}' > ${mut}-2.covThreshold${t}.neg.bed

    # Count number of reads in possible pauses for both replicates
    reads1=`cat ${mut}-1.covThreshold${t}.???.bed | awk -F "\t" '{sum += $7} END {print sum}'`
    reads2=`cat ${mut}-2.covThreshold${t}.???.bed | awk -F "\t" '{sum += $7} END {print sum}'`

    # Count number of reads in pauses for both replicates
    pauses1=`cat ${mut}_PASS_pauseScore_cov${t}.txt | awk -F "\t" '{sum += $4} END {print sum}'`
    pauses2=`cat ${mut}_PASS_pauseScore_cov${t}.txt | awk -F "\t" '{sum += $5} END {print sum}'` 

    echo -e "${mut}\t${reads1}\t${reads2}\t${pauses1}\t${pauses2}" >> pauseStrength_cov${t}.txt

done < $1

./Scripts/plotPauseStrength.R $t
