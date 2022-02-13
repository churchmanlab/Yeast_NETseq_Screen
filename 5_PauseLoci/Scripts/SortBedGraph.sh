#!/bin/bash

#SBATCH -c 1
#SBATCH -t 0-00:30
#SBATCH -p short
#SBATCH --mem=1G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=youremail@gmail.com


### Use
# samples="cac1-1 cac1-2 cac2-1 cac2-2 cac3-1 cac3-2"
# # # samples="wt-1 wt-2 wt-3 wt-4 dst1-1 dst1-2 swr1-1 swr1-2 bre1-1 bre1-2 bur2-1 bur2-2 bye1-1 bye1-2 cbc1-1 cbc1-2 ccr4-1 ccr4-2 cdc39-1 cdc39-2 cdc73-1 cdc73-2 chd1-1 chd1-2 ctr9-1 ctr9-2 dhh1-1 eaf1-1 eaf1-2 elf1-1 elf1-2 gcn5-1 gcn5-2 hda1-1 hda1-2 hpc2-1 hpc2-2 htz1-1 htz1-2 ino80-1 ino80-2 isw1-1 isw1-2 isw2-1 isw2-2 leo1-1 leo1-2 nap1-1 nap1-2 npl3-1 npl3-2 paf1-1 paf1-2 rad6-1 rad6-2 rco1-2 rpb4-1 rpb4-2 rph1-1 rph1-2 rsc30-1 rsc30-2 rtf1-1 rtf1-2 rtt103-1 rtt103-2 set2-1 set2-2 set3-1 set3-2 spt4-1 spt4-2 ubp8-1 ubp8-2 vps15-1 vps15-2 vps34-1 vps34-2"
# for s in $samples
# do
# sbatch -o logs/SortBG_${sample}.log -e logs/SortBG_${sample}.err ./Scripts/SortBedGraph.sh ${s}_noDups_stt_pos.bedGraph ${s}_noDups_stt_neg.bedGraph
# sbatch -o logs/SortBG_${sample}.log -e logs/SortBG_${sample}.err ./Scripts/SortBedGraph.sh ${s}_stt_pos.bedGraph ${s}_stt_neg.bedGraph
# done

###

Pos=$1
Neg=$2

# Get bedgraph files
posBG=$Pos
negBG=$Neg

# Remove mtDNA mappers, convert roman to arabic numerals, sort, and convert back
sed '/^chrM/d' $posBG | sed -e $'s/chrI\t/chr1\t/' -e $'s/chrII\t/chr2\t/' -e $'s/chrIII\t/chr3\t/' -e $'s/chrIV\t/chr4\t/' -e $'s/chrV\t/chr5\t/' -e $'s/chrVI\t/chr6\t/' -e $'s/chrVII\t/chr7\t/' -e $'s/chrVIII\t/chr8\t/' -e $'s/chrIX\t/chr9\t/' -e $'s/chrX\t/chr10\t/' -e $'s/chrXI\t/chr11\t/' -e $'s/chrXII\t/chr12\t/' -e $'s/chrXIII\t/chr13\t/' -e $'s/chrXIV\t/chr14\t/' -e $'s/chrXV\t/chr15\t/' -e $'s/chrXVI\t/chr16\t/' | sort -nk1.4,1 -k2,2n | sed -e $'s/chr1\t/chrI\t/' -e $'s/chr2\t/chrII\t/' -e $'s/chr3\t/chrIII\t/' -e $'s/chr4\t/chrIV\t/' -e $'s/chr5\t/chrV\t/' -e $'s/chr6\t/chrVI\t/' -e $'s/chr7\t/chrVII\t/' -e $'s/chr8\t/chrVIII\t/' -e $'s/chr9\t/chrIX\t/' -e $'s/chr10\t/chrX\t/' -e $'s/chr11\t/chrXI\t/' -e $'s/chr12\t/chrXII\t/' -e $'s/chr13\t/chrXIII\t/' -e $'s/chr14\t/chrXIV\t/' -e $'s/chr15\t/chrXV\t/' -e $'s/chr16\t/chrXVI\t/' > ${posBG/.bed/_sort.bed}

sed '/^chrM/d' $negBG | sed -e $'s/chrI\t/chr1\t/' -e $'s/chrII\t/chr2\t/' -e $'s/chrIII\t/chr3\t/' -e $'s/chrIV\t/chr4\t/' -e $'s/chrV\t/chr5\t/' -e $'s/chrVI\t/chr6\t/' -e $'s/chrVII\t/chr7\t/' -e $'s/chrVIII\t/chr8\t/' -e $'s/chrIX\t/chr9\t/' -e $'s/chrX\t/chr10\t/' -e $'s/chrXI\t/chr11\t/' -e $'s/chrXII\t/chr12\t/' -e $'s/chrXIII\t/chr13\t/' -e $'s/chrXIV\t/chr14\t/' -e $'s/chrXV\t/chr15\t/' -e $'s/chrXVI\t/chr16\t/' | sort -nk1.4,1 -k2,2n | sed -e $'s/chr1\t/chrI\t/' -e $'s/chr2\t/chrII\t/' -e $'s/chr3\t/chrIII\t/' -e $'s/chr4\t/chrIV\t/' -e $'s/chr5\t/chrV\t/' -e $'s/chr6\t/chrVI\t/' -e $'s/chr7\t/chrVII\t/' -e $'s/chr8\t/chrVIII\t/' -e $'s/chr9\t/chrIX\t/' -e $'s/chr10\t/chrX\t/' -e $'s/chr11\t/chrXI\t/' -e $'s/chr12\t/chrXII\t/' -e $'s/chr13\t/chrXIII\t/' -e $'s/chr14\t/chrXIV\t/' -e $'s/chr15\t/chrXV\t/' -e $'s/chr16\t/chrXVI\t/' > ${negBG/.bed/_sort.bed}

# Make executable for later step
chmod u+rwx ${negBG/.bed/_sort.bed}
chmod u+rwx ${posBG/.bed/_sort.bed}


