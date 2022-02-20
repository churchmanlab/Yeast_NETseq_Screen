#!/bin/bash

#SBATCH -c 1
#SBATCH -t 0-12:00
#SBATCH -p short
#SBATCH --mem=30G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=youremail@gmail.com


# Use: 
# samples="bre1 bur2 bye1 cbc1 ccr4 cdc39 cdc73 chd1 ctr9 dst1 eaf1 elf1 hda1 hpc2 htz1 ino80 isw1 isw2 leo1 nap1 npl3 paf1 rad6 rco1 rpb4 rph1 rsc30 rtf1 rtt103 set2 set3 spt4 swr1 ubp8 vps15 vps34 wt"
# types="Real Shuffled"
# for type in $types
# do
# for sample in $samples
# do
# sbatch -o logs/6_calcFeatures${type}_${sample}.log -e logs/6_calcFeatures${type}_${sample}.err ./Scripts/run_calcFeatures_MC.sh $sample $type
# done
# done

# After all jobs complete:
# for sample in $samples
# do
# cat ${sample}.IDRrep.PAUSESwithFEATURES_*.bed > ${sample}.IDRrep.PAUSESwithFEATURES.bed
# done





mut=$1
type=$2

./Scripts/calcFeatures_MC.R $mut $type
