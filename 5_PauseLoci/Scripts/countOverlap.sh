#!/bin/bash

#SBATCH -c 1
#SBATCH -t 0-01:00
#SBATCH -p short
#SBATCH --mem=10G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=youremail@gmail.com
###
# Use: 
# t=2
# sbatch -o logs/5_countOverlap.log -e logs/5_countOverlap.err ./Scripts/countOverlap.sh deletionStrainsCovFiltered.txt $t
# 
t=$2

while read mut
do

    echo $mut

    start=`awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, "pause"NR, ".", $4}' ${mut}.IDRrepPausesAll_cov${t}.bed | intersectBed -a ../0_Annotations/Genes/txptsStart15p.bed -b - -wb -s | wc -l`
    mid=`awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, "pause"NR, ".", $4}' ${mut}.IDRrepPausesAll_cov${t}.bed | intersectBed -a ../0_Annotations/Genes/txptsMid.bed -b - -wb -s | wc -l`
    end=`awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, "pause"NR, ".", $4}' ${mut}.IDRrepPausesAll_cov${t}.bed | intersectBed -a ../0_Annotations/Genes/txptsEnd85p.bed -b - -wb -s | wc -l`

    total=$(echo $start + $mid + $end  | bc)

    echo -e "$mut\t$start\t$mid\t$end\t$total" >> geneOverlap.txt

done < $1

# Shuffled control
start=`awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, "pause"NR, ".", $4}' ALL.IDRrep.SHUFFLE.bed | intersectBed -a ../0_Annotations/Genes/txptsStart15p.bed -b - -wb -s | wc -l`
mid=`awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, "pause"NR, ".", $4}' ALL.IDRrep.SHUFFLE.bed | intersectBed -a ../0_Annotations/Genes/txptsMid.bed -b - -wb -s | wc -l`
end=`awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, "pause"NR, ".", $4}' ALL.IDRrep.SHUFFLE.bed | intersectBed -a ../0_Annotations/Genes/txptsEnd85p.bed -b - -wb -s | wc -l`

total=$(echo $start + $mid + $end  | bc)

echo -e "zzzCNTL\t$start\t$mid\t$end\t$total" >> geneOverlap.txt
