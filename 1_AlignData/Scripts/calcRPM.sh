#!/bin/bash

# Normalize samples by reads per million

# Written by: K. Lachance
# Date: November 26, 2018

# Use: ./combineReps.sh mut

# Input parameters
mut=$1

# Make output directory
mkdir -p ./RPMnormFiles

# Calculate coverage on each strand
sumPos=`cat ./CoverageFiles/${mut}_stt_pos.bedGraph | grep -v track | awk -F "\t" '{sum += ($3-$2)*$4} END {printf "%.0f", sum}'`
sumNeg=`cat ./CoverageFiles/${mut}_stt_neg.bedGraph | grep -v track | awk -F "\t" '{sum += ($3-$2)*$4} END {printf "%.0f", sum}'`

# Calculate total coverage
cov=$(($sumPos + $sumNeg))
echo "${mut}: ${cov}"

# Calculate normalization factor
norm=`echo $cov/1000000 | bc -l`

# Create finished bed files (RPM normalized)
awk -F "\t" -v norm=$norm '{print $1, $2, $3, $4/(norm)}' ./CoverageFiles/${mut}*_stt_pos.bedGraph | sed 's/ /\t/g'> ./RPMnormFiles/${mut}.pos.bed
awk -F "\t" -v norm=$norm '{print $1, $2, $3, $4/(norm)}' ./CoverageFiles/${mut}*_stt_neg.bedGraph | sed 's/ /\t/g' > ./RPMnormFiles/${mut}.neg.bed

