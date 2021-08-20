#!/bin/bash

#SBATCH -c 1
#SBATCH -t 0-01:00
#SBATCH -p short
#SBATCH --mem=10G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=mtcouvi@gmail.com


### Use: 
# samples="wt dst1 swr1"
samples="bre1 bur2 bye1 cbc1 ccr4 cdc39 cdc73 chd1 ctr9 dst1 eaf1 elf1 hda1 hpc2 htz1 ino80 isw1 isw2 leo1 nap1 npl3 paf1 rad6 rco1 rpb4 rph1 rsc30 rtf1 rtt103 set2 set3 spt4 swr1 ubp8 vps15 vps34 wt"
# samples="spt4 npl3 htz1"
t=2
for sample in $samples
do
sbatch -o logs/5_shufflePausesfasta_${sample}.log -e logs/5_shufflePausesfasta_${sample}.err ./Scripts/shufflePauses_MC.sh ${sample} ${t}
done

mut=$1
t=$2

# Genome file for chr ordering
genome=../0_Annotations/Scerevisiae_R64/SacCer3_R64_genomic.fasta.fai

# For real pauses
# cat ${mut}_PASS_pauseScore_cov${t}.txt | awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, $6}' > ${mut}.IDRrep.PAUSES.bed 
# cat ${mut}_PASS_pauseScore_cov${t}.txt | awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2-10, $3+10, $6}' > ${mut}.IDRrep.pauseWindows.bed
# bedtools getfasta -fi ../0_Annotations/Scerevisiae_R64/SacCer3_R64_genomic.fasta -bed ${mut}.IDRrep.pauseWindows.bed -s > ${mut}.IDRrep.pauseWindows.fasta
# 
# ## #To shuffle pauses
# ## #cat ${mut}.covThreshold${t}.pos.bed ${mut}.covThreshold${t}.neg.bed | sortBed -g $genome -i > ${mut}.covThreshold${t}.bed
# ## cat ${mut}.covThreshold${t}.bed | grep "+" > ${mut}.covThreshold${t}.pos.bed
# ## cat ${mut}.covThreshold${t}.bed | grep -v "+" > ${mut}.covThreshold${t}.neg.bed
# cat ${mut}.IDRrep.PAUSES.bed | grep "+" | sortBed -i - > ${mut}.IDRrep.PAUSES.pos.bed 
# cat ${mut}.IDRrep.PAUSES.bed | grep -v "+" | sortBed -i - > ${mut}.IDRrep.PAUSES.neg.bed
# 
# mergeBed -i ${mut}.covThreshold${t}.pos.bed | intersectBed -a - -b ${mut}.IDRrep.PAUSES.pos.bed -c > ${mut}.pausesByGene.pos.bed
# mergeBed -i ${mut}.covThreshold${t}.neg.bed | intersectBed -a - -b ${mut}.IDRrep.PAUSES.neg.bed -c > ${mut}.pausesByGene.neg.bed
# 
# ./Scripts/shufflePauses.R ${mut} pos
# ./Scripts/shufflePauses.R ${mut} neg
# 
# cat ${mut}.IDRrep.SHUFFLE.pos.bed ${mut}.IDRrep.SHUFFLE.neg.bed | sortBed -g $genome -i > ${mut}.IDRrep.SHUFFLE.bed

cat ${mut}.IDRrep.SHUFFLE.pos.bed ${mut}.IDRrep.SHUFFLE.neg.bed | awk -F "\t" 'BEGIN {OFS="\t"} {print $1, $2-10, $3+10, $6}' > ${mut}.IDRrep.pauseWindowsSHUFFLE.bed
bedtools getfasta -fi ../0_Annotations/Scerevisiae_R64/SacCer3_R64_genomic.fasta -bed ${mut}.IDRrep.pauseWindowsSHUFFLE.bed -s > ${mut}.IDRrep.pauseWindowsSHUFFLE.fasta


# Organize files
#mkdir -p ${mut}
#mv ${mut}*bed ${mut}*fasta ${mut}

# if error, check ${mut}.IDRrep.SHUFFLE.pos.bed and ${mut}.IDRrep.SHUFFLE.neg.bed for scientific notation