#!/bin/bash

#SBATCH -c 1
#SBATCH -t 0-12:00
#SBATCH -p short
#SBATCH --mem=50G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=mtcouvi@gmail.com

### Script to identify pause locations in genes using NET-seq data (as defined by Churchman & Weissman, 2011)

# Written by: K. Lachance
# Date: November 1, 2018
# Updated by: M. Couvillion 5/2021
# - To use updated findPauses_NB.sh (findPauses_NB_MC.sh) that removes gene "buffer" at beginning and end that excludes those regions from pause calling. Instead make window up/down size change if close to ends, to keep total window size constant

# Use: 
# samples="bre1 bur2 bye1 cbc1 ccr4 cdc39 cdc73 chd1 ctr9 dhh1 dst1 eaf1 elf1 gcn5 hda1 hpc2 htz1 ino80 isw1 isw2 leo1 nap1 npl3 paf1 rad6 rco1 rpb4 rph1 rsc30 rtf1 rtt103 set2 set3 spt4 swr1 ubp8 vps15 vps34 wt"
# t=2
# strands="pos neg" 
# for strand in $strands
# do
# for sample in $samples
# do
# sbatch -o logs/5_findPauses_cov${t}_${sample}-1_${strand}.log -e logs/5_findPauses_cov${t}_${sample}-1_${strand}.err ./Scripts/run_findPausesRepCov.sh $sample ${sample}-1 $t $strand
# 
# sbatch -o logs/5_findPauses_cov${t}_${sample}-2_${strand}.log -e logs/5_findPauses_cov${t}_${sample}-2_${strand}.err ./Scripts/run_findPausesRepCov.sh $sample ${sample}-2 $t $strand
# done
# done


module load bedops/2.4.30
module load bedtools/2.27.1


# Mutant of interest
mut=$1
rep=$2
# Txpt coverage threshold
t=$3
# Strand
strand=$4

# Echo deletion strain
echo $mut

# standard deviation over mean to call as pause
SD="3" 
# Window for NB distribution (actually half of window size - nts on each side of potential pause)
w=100

# Get bedgraph files
BgPath="../1_AlignData/"
suffix=stt_${strand}_sort.bedGraph
Signal=${BgPath}${rep}_${suffix} 


# Redefine gene files
Genes=${mut}.covThreshold${t}.${strand}.bed

# Counter for progress
pos=`wc -l $Genes | awk '{print $1}'`
i=1

echo ${rep} >> logs/5_findPauses_cov${t}_${rep}_${strand}.err
# For each gene (pos)
while read gene
do
	echo $gene >> logs/5_findPauses_cov${t}_${rep}_${strand}.err
	echo $gene
	
    # Progress tracker
    if (( $i % 100 == 0 )) 
    then
	echo "${i} / ${pos} genes done!"
    fi

    # Store gene as one-line to temporary file
    echo $gene > ${rep}_cov${t}_${strand}.tmpGene.bed

    # Overlap gene with NET-seq signal, with one row per basepair
    bedops --chop ${rep}_cov${t}_${strand}.tmpGene.bed | intersectBed -a $Signal -b - > ${rep}_cov${t}_${strand}.tmpCov.bed

    # Call R script to find pauses
#     rm ${mut}_cov${t}.pausesNB_T${SD}.w${w}.${strand}.bed
    
    ./Scripts/findPauses_NB_MC.R ${rep}_cov${t}_${strand}.tmpCov.bed ${rep} ${t} ${SD} ${w} ${strand}

    #  Increment tracker
    i=$((i+1))

done < $Genes
echo "Done!"

