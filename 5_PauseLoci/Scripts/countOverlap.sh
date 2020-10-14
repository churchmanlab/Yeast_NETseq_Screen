#!/bin/bash

while read mut
do

    echo $mut

    start=`awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, "pause"NR, ".", $4}' ${mut}.IDRrepPausesAll.bed | intersectBed -a ../0_Annotations/Genes/geneStart15p.bed -b - -wb -s | wc -l`
    mid=`awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, "pause"NR, ".", $4}' ${mut}.IDRrepPausesAll.bed | intersectBed -a ../0_Annotations/Genes/geneMid.bed -b - -wb -s | wc -l`
    end=`awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, "pause"NR, ".", $4}' ${mut}.IDRrepPausesAll.bed | intersectBed -a ../0_Annotations/Genes/geneEnd85p.bed -b - -wb -s | wc -l`

    total=$(echo $start + $mid + $end  | bc)

    echo -e "$mut\t$start\t$mid\t$end\t$total" >> geneOverlap.txt

done < $1

# Shuffled control
start=`awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, "pause"NR, ".", $4}' ALL.IDRrep.SHUFFLE.bed | intersectBed -a ../0_Annotations/Genes/geneStart15p.bed -b - -wb -s | wc -l`
mid=`awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, "pause"NR, ".", $4}' ALL.IDRrep.SHUFFLE.bed | intersectBed -a geneMid.bed -b - -wb -s | wc -l`
end=`awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, "pause"NR, ".", $4}' ALL.IDRrep.SHUFFLE.bed | intersectBed -a ../0_Annotations/Genes/geneEnd85p.bed -b - -wb -s | wc -l`

total=$(echo $start + $mid + $end  | bc)

echo -e "zzzCNTL\t$start\t$mid\t$end\t$total" >> geneOverlap.txt
