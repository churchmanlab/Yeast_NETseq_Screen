#!/bin/bash

#SBATCH -c 1
#SBATCH -t 0-03:00
#SBATCH -p short
#SBATCH --mem=5G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=mtcouvi@gmail.com

### Script to choose genes passing coverage threshold in both replicates for pause calling

# Written by: K. Lachance
# Date: November 1, 2018
# Updated by: M. Couvillion 5/2021
# - The original run_findPauses.sh script was split into this one and run_findPausesRepCov.sh
# - Use transcript bounds instead of CDS
# - Filter out ncRNA genes

# Use: 
# samples="bre1 bur2 bye1 cbc1 ccr4 cdc39 cdc73 chd1 ctr9 dhh1 dst1 eaf1 elf1 gcn5 hda1 hpc2 htz1 ino80 isw1 isw2 leo1 nap1 npl3 paf1 rad6 rco1 rpb4 rph1 rsc30 rtf1 rtt103 set2 set3 spt4 swr1 ubp8 vps15 vps34 wt"
# mkdir logs
# t=2
# strands="pos neg" 
# for strand in $strands
# do
# for sample in $samples
# do
# sbatch -o logs/5_findGenes_cov${t}_${sample}_${strand}.log -e logs/5_findGenes_cov${t}_${sample}_${strand}.err ./Scripts/run_findGenesOverCovThreshold_MC.sh ${sample} $t $strand
# done
# done


module load bedops/2.4.30
module load bedtools/2.27.1


# Mutant of interest
mut=$1
# Txpt coverage threshold
t=$2
# Strand
strand=$3

# Echo deletion strain
echo $mut

BgPath="../1_AlignData/"

suffix=stt_${strand}_sort.bedGraph



# Get gene files
Genes=../0_Annotations/Genes/txpts_all_${strand}.bed

reps="1 2"
for rep in $reps
do

# Get bedgraph files
Signal=${BgPath}${mut}-${rep}_${suffix} 

# Overlap and isolate genes that are above coverage threshold
rm ${mut}-${rep}.covThreshold${t}.${strand}.bed

sed '1d' $Signal | awk -F "\t" 'BEGIN {OFS="\t"} ($4 > 0) {print $1, $2, $3, $4*($3-$2)}' | mapBed -a $Genes -b - -c 4 -o sum | awk -F "\t" 'BEGIN {OFS="\t"} ($7 != ".") {print $0, $7/($3-$2)}' | awk -F "\t" -v t=$t 'BEGIN {OFS="\t"} ($8 > t) {print $1, $2, $3, $4, $5, $6}' > ${mut}-${rep}.covThreshold${t}.${strand}.bed

done

# Only take overlapping genes from two replicates
grep -f ${mut}-1.covThreshold${t}.${strand}.bed ${mut}-2.covThreshold${t}.${strand}.bed > ${mut}.covThreshold${t}.${strand}.bed

