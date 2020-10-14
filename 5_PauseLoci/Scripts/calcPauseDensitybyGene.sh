#!/bin/bash

while read mut
do

    echo $mut

    # Format pause file properly
    cat ${mut}_PASS_pauseScore.txt | awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, $6}' > ${mut}.IDRrepPausesAll.bed
    
    # Calculate pause density
    awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, "pause"NR, ".", $4}' ${mut}.IDRrepPausesAll.bed | intersectBed -a ${mut}.covThreshold5.bed  -b - -c -s | awk -F "\t" -v mut=${mut} 'BEGIN {OFS="\t"} {print $0, ($7*1000)/($3-$2), mut}' > ${mut}.pauseDensity.txt

done < $1

# Combine all files together
cat *.pauseDensity.txt > ALL.pauseDensity.txt

# Plot pause density
./Scripts/plot_pauseDensity.R
