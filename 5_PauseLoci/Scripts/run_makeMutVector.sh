#!/bin/bash

# Combine all pause locations into single file of all uniq pauses 
cat *_PASS_pauseScore.txt | cut -f 1,2,3,6 | sort | uniq > AllPauseLoci.bed

while read mut
do
    echo $mut
    intersectBed -a AllPauseLoci.bed -b ${mut}.covThreshold5.* > ${mut}.possiblePauses.bed
    ./Scripts/makeMutVector.R AllPauseLoci.bed $mut

done < $1
