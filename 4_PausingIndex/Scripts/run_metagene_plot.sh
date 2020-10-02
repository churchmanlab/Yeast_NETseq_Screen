#!/bin/bash

mkdir -p readyToPlot

while read mut
do

    echo $mut

    awk '($4 != 0) ' ${mut}.TSS_OverlapMatrix.350.pos.txt > tmp.pos.txt
    awk '($4 != 0) ' ${mut}.TSS_OverlapMatrix.350.neg.txt> tmp.neg.txt
    ./Scripts/metagene_plot.R tmp.pos.txt tmp.neg.txt ${mut} TSS

    awk '($4 != 0) ' ${mut}.pA_OverlapMatrix.350.pos.txt > tmp.pos.txt
    awk '($4 != 0) ' ${mut}.pA_OverlapMatrix.350.neg.txt> tmp.neg.txt
    ./Scripts/metagene_plot.R tmp.pos.txt tmp.neg.txt ${mut} pA

    awk '($4 != 0) ' ${mut}.anti_OverlapMatrix.350.pos.txt > tmp.pos.txt
    awk '($4 != 0) ' ${mut}.anti_OverlapMatrix.350.neg.txt> tmp.neg.txt
    ./Scripts/metagene_plot.R tmp.pos.txt tmp.neg.txt ${mut} Anti

    awk '($4 != 0) ' ${mut}.5pSS_OverlapMatrix.25.pos.txt > tmp.pos.txt
    awk '($4 != 0) ' ${mut}.5pSS_OverlapMatrix.25.neg.txt> tmp.neg.txt
    ./Scripts/metagene_plot.R tmp.pos.txt tmp.neg.txt ${mut} 5pSS

    awk '($4 != 0) ' ${mut}.3pSS_OverlapMatrix.25.pos.txt > tmp.pos.txt
    awk '($4 != 0) ' ${mut}.3pSS_OverlapMatrix.25.neg.txt> tmp.neg.txt
    ./Scripts/metagene_plot.R tmp.pos.txt tmp.neg.txt ${mut} 3pSS

    rm -fr tmp.pos.txt tmp.neg.txt

#done < deletionStrains.txt
done < tmp
