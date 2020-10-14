#!/bin/bash

# Combine and normalize highly correlated replicates

# Written by: K. Lachance
# Date: November 26, 2018

# Use: ./combineReps.sh mut

# Input parameters
mut=$1

# Calculate number of replicates
numReps=`ls ./CoverageFiles/${mut}*_stt_pos.bedGraph | wc -l`
echo "There are ${numReps} replicates for ${mut}"

# Make output directory
mkdir -p ./CombinedReplicates

# Combine unnormalized replicates (numReps = 1)
if [ $numReps -eq 1 ]
then
    cp ./CoverageFiles/${mut}*_stt_pos.bedGraph ./CombinedReplicates/${mut}.noNorm.pos.bed
    cp ./CoverageFiles/${mut}*_stt_neg.bedGraph ./CombinedReplicates/${mut}.noNorm.neg.bed    
fi

# Combine unnormalized replicates (numReps = 2)
if [ $numReps -eq 2 ]
then
    bedtools unionbedg -i ./CoverageFiles/${mut}*_stt_pos.bedGraph | awk -v numReps=$numReps -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, ($4+$5)}' > ./CombinedReplicates/${mut}.noNorm.pos.bed
    bedtools unionbedg -i ./CoverageFiles/${mut}*_stt_neg.bedGraph | awk -v numReps=$numReps -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, ($4+$5)}' > ./CombinedReplicates/${mut}.noNorm.neg.bed
fi

# Combine unnormalized replicates (numReps = 3)
if [ $numReps -eq 3 ]
then
    bedtools unionbedg -i ./CoverageFiles/${mut}*_stt_pos.bedGraph | awk -v numReps=$numReps -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, ($4+$5+$6)}' > ./CombinedReplicates/${mut}.noNorm.pos.bed
    bedtools unionbedg -i ./CoverageFiles/${mut}*_stt_neg.bedGraph | awk -v numReps=$numReps -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, ($4+$5+$6)}' > ./CombinedReplicates/${mut}.noNorm.neg.bed
fi

# Combine unnormalized replicates (numReps = 4)
if [ $numReps -eq 4 ]
then
    bedtools unionbedg -i ./CoverageFiles/${mut}*_stt_pos.bedGraph | awk -v numReps=$numReps -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, ($4+$5+$6+$7)}' > ./CombinedReplicates/${mut}.noNorm.pos.bed
    bedtools unionbedg -i ./CoverageFiles/${mut}*_stt_neg.bedGraph | awk -v numReps=$numReps -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, ($4+$5+$6+$7)}' > ./CombinedReplicates/${mut}.noNorm.neg.bed
fi

# Combine unnormalized replicates (numReps = 5)
if [ $numReps -eq 5 ]
then
    bedtools unionbedg -i ./CoverageFiles/${mut}*_stt_pos.bedGraph | awk -v numReps=$numReps -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, ($4+$5+$6+$7+$8)}' > ./CombinedReplicates/${mut}.noNorm.pos.bed
    bedtools unionbedg -i ./CoverageFiles/${mut}*_stt_neg.bedGraph | awk -v numReps=$numReps -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, ($4+$5+$6+$7+$8)}' > ./CombinedReplicates/${mut}.noNorm.neg.bed
fi

# Combine unnormalized replicates (numReps = 6)
if [ $numReps -eq 6 ]
then
    bedtools unionbedg -i ./CoverageFiles/${mut}*_stt_pos.bedGraph | awk -v numReps=$numReps -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, ($4+$5+$6+$7+$8+$9)}' > ./CombinedReplicates/${mut}.noNorm.pos.bed
    bedtools unionbedg -i ./CoverageFiles/${mut}*_stt_neg.bedGraph | awk -v numReps=$numReps -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, ($4+$5+$6+$7+$8+$9)}' > ./CombinedReplicates/${mut}.noNorm.neg.bed
fi

# Calculate coverage on each strand
sumPos=`awk -F "\t" '{sum += ($3-$2)*$4} END {printf "%.0f", sum}' ./CombinedReplicates/${mut}.noNorm.pos.bed`
sumNeg=`awk -F "\t" '{sum += ($3-$2)*$4} END {printf "%.0f", sum}' ./CombinedReplicates/${mut}.noNorm.neg.bed`

# Calculate total coverage
cov=$(($sumPos + $sumNeg))
echo "${mut}: ${cov}"

# Calculate normalization factor
norm=`echo $cov/1000000 | bc -l`

# Create finished bed files (RPM and numRep normalized)
awk -F "\t" -v norm=$norm -v numReps=$numReps '{print $1, $2, $3, $4/(norm*numReps)}' ./CombinedReplicates/${mut}.noNorm.pos.bed | sed 's/ /\t/g'> ./CombinedReplicates/${mut}.pos.bed
awk -F "\t" -v norm=$norm -v numReps=$numReps '{print $1, $2, $3, $4/(norm*numReps)}' ./CombinedReplicates/${mut}.noNorm.neg.bed | sed 's/ /\t/g' > ./CombinedReplicates/${mut}.neg.bed

