#!/bin/bash

#SBATCH -c 1
#SBATCH -t 0-01:00
#SBATCH -p short
#SBATCH --mem=10G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<youremail@address.com>
###
# Use: 
# covThreshold=2
# sbatch -o logs/5_calcPDbyGene_${covThreshold}.log -e logs/5_calcPDbyGene_${covThreshold}.err ./Scripts/calcPauseDensitybyGene.sh deletionStrains.txt $covThreshold


# Txpt coverage threshold
t=$2

while read mut
do

    echo $mut
    
    # Combine covThreshold beds (must be present in both files)
	grep -f ${mut}-1.covThreshold${t}.pos.bed ${mut}-2.covThreshold${t}.pos.bed >  ${mut}.covThreshold${t}.pos.bed
	grep -f ${mut}-1.covThreshold${t}.neg.bed ${mut}-2.covThreshold${t}.neg.bed >  ${mut}.covThreshold${t}.neg.bed
	cat ${mut}.covThreshold${t}.pos.bed ${mut}.covThreshold${t}.neg.bed > ${mut}.covThreshold${t}.bed
	
    # Format pause file properly
    cat ${mut}_PASS_pauseScore_cov${t}.txt | awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, $6}' > ${mut}.IDRrepPausesAll_cov${t}.bed
    
    # Calculate pause density
    awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, "pause"NR, ".", $4}' ${mut}.IDRrepPausesAll_cov${t}.bed | intersectBed -a ${mut}.covThreshold${t}.bed -b - -c -s | awk -F "\t" -v mut=${mut} 'BEGIN {OFS="\t"} {print $0, ($7*1000)/($3-$2), mut}' > ${mut}_cov${t}.pauseDensity.txt

done < $1

# Combine all files together
cat *_cov${t}.pauseDensity.txt > ALL_cov${t}.pauseDensity.txt

# Plot pause density
./Scripts/plot_pauseDensity.R $t
