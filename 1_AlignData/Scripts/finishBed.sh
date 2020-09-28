#!/bin/bash

# Written by: K. Lachance

# Date: November 26, 2018

# Script to calculate total # of reads and convert .stt.bed files to final, sorted bedgraph format

# Input
strain=$1 # Fastq file name (e.g. BRE1_ES_216)
mut=$2 # name of mutant replicate (e.g. bre1-1)

# Calculate coverage on each strand
sumPos=`awk -F "\t" '{sum += ($3-$2)*$4} END {print sum}' ${strain}_stt_pos.bedGraph`
sumNeg=`awk -F "\t" '{sum += ($3-$2)*$4} END {print sum}' ${strain}_stt_neg.bedGraph`

# Calculate total coverage
cov=$(($sumPos + $sumNeg))
echo "${strain}: ${cov}"

# Calculate normalization factor
norm=`echo $cov/1000000 | bc -l`

# Create finished bed files (RPM normalized)
sortBed -i ${strain}_stt_pos.bedGraph | awk -F "\t" -v norm=$norm '{print $1, $2, $3, $4/norm}' | sed 's/ /\t/g'> ../DATA/RPMnormBedgraphFiles/${mut}.pos.bed
sortBed -i ${strain}_stt_neg.bedGraph | awk -F "\t" -v norm=$norm '{print $1, $2, $3, $4/norm}' | sed 's/ /\t/g' > ../DATA/RPMnormBedgraphFiles/${mut}.neg.bed

# Move other coverage files
mv ${strain}*.bedGraph ./CompletedCoverageFiles/
sortBed -i ./CompletedCoverageFiles/${strain}_stt_pos.bedGraph > ../DATA/NonormBedgraphFiles/${mut}.pos.bed
sortBed -i ./CompletedCoverageFiles/${strain}_stt_neg.bedGraph > ../DATA/NonormBedgraphFiles/${mut}.neg.bed
