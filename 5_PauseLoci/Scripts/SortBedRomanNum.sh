#!/bin/bash

#SBATCH -c 1
#SBATCH -t 0-00:30
#SBATCH -p short
#SBATCH --mem=1G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<youremail@address.com>


### Use: 
# sbatch ./Scripts/SortBedRomanNum.sh rep

rep=$1 


# Pos='txpts2013_all_pos.bed'
# Neg='txpts2013_all_neg.bed'
Pos=../1_AlignData/${rep}_stt_pos.bedGraph
Neg=../1_AlignData/${rep}_stt_neg.bedGraph
# Pos=${rep}.covThreshold5.pos.bed
# Neg=${rep}.covThreshold5.neg.bed



# Remove mtDNA mappers, convert roman to arabic numerals, sort, and convert back
sed '/^chrM/d' $Pos | sed -e $'s/chrI\t/chr1\t/' -e $'s/chrII\t/chr2\t/' -e $'s/chrIII\t/chr3\t/' -e $'s/chrIV\t/chr4\t/' -e $'s/chrV\t/chr5\t/' -e $'s/chrVI\t/chr6\t/' -e $'s/chrVII\t/chr7\t/' -e $'s/chrVIII\t/chr8\t/' -e $'s/chrIX\t/chr9\t/' -e $'s/chrX\t/chr10\t/' -e $'s/chrXI\t/chr11\t/' -e $'s/chrXII\t/chr12\t/' -e $'s/chrXIII\t/chr13\t/' -e $'s/chrXIV\t/chr14\t/' -e $'s/chrXV\t/chr15\t/' -e $'s/chrXVI\t/chr16\t/' | sort -nk1.4,1 -k2,2n | sed -e $'s/chr1\t/chrI\t/' -e $'s/chr2\t/chrII\t/' -e $'s/chr3\t/chrIII\t/' -e $'s/chr4\t/chrIV\t/' -e $'s/chr5\t/chrV\t/' -e $'s/chr6\t/chrVI\t/' -e $'s/chr7\t/chrVII\t/' -e $'s/chr8\t/chrVIII\t/' -e $'s/chr9\t/chrIX\t/' -e $'s/chr10\t/chrX\t/' -e $'s/chr11\t/chrXI\t/' -e $'s/chr12\t/chrXII\t/' -e $'s/chr13\t/chrXIII\t/' -e $'s/chr14\t/chrXIV\t/' -e $'s/chr15\t/chrXV\t/' -e $'s/chr16\t/chrXVI\t/' > ${Pos/.bed/_sort.bed}

sed '/^chrM/d' $Neg | sed -e $'s/chrI\t/chr1\t/' -e $'s/chrII\t/chr2\t/' -e $'s/chrIII\t/chr3\t/' -e $'s/chrIV\t/chr4\t/' -e $'s/chrV\t/chr5\t/' -e $'s/chrVI\t/chr6\t/' -e $'s/chrVII\t/chr7\t/' -e $'s/chrVIII\t/chr8\t/' -e $'s/chrIX\t/chr9\t/' -e $'s/chrX\t/chr10\t/' -e $'s/chrXI\t/chr11\t/' -e $'s/chrXII\t/chr12\t/' -e $'s/chrXIII\t/chr13\t/' -e $'s/chrXIV\t/chr14\t/' -e $'s/chrXV\t/chr15\t/' -e $'s/chrXVI\t/chr16\t/' | sort -nk1.4,1 -k2,2n | sed -e $'s/chr1\t/chrI\t/' -e $'s/chr2\t/chrII\t/' -e $'s/chr3\t/chrIII\t/' -e $'s/chr4\t/chrIV\t/' -e $'s/chr5\t/chrV\t/' -e $'s/chr6\t/chrVI\t/' -e $'s/chr7\t/chrVII\t/' -e $'s/chr8\t/chrVIII\t/' -e $'s/chr9\t/chrIX\t/' -e $'s/chr10\t/chrX\t/' -e $'s/chr11\t/chrXI\t/' -e $'s/chr12\t/chrXII\t/' -e $'s/chr13\t/chrXIII\t/' -e $'s/chr14\t/chrXIV\t/' -e $'s/chr15\t/chrXV\t/' -e $'s/chr16\t/chrXVI\t/' > ${Neg/.bed/_sort.bed}

# Make executable for later step
# chmod u+rwx ${Pos/.bed/_sort.bed}
# chmod u+rwx ${Neg/.bed/_sort.bed}


