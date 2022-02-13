#!/bin/bash

#SBATCH -c 1
#SBATCH -t 0-01:00
#SBATCH -p short
#SBATCH --mem=1G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=mtcouvi@gmail.com

### Get read totals for each strain from bedGraph files

# Use: 
# echo -e 'wt\ndst1\nswr1' > tmpStrains.txt
# sbatch -o logs/5_GetReadTotals.log -e logs/5_GetReadTotals.err ./Scripts/GetReadTotals_MC.sh

rm readTotals.txt
rm readTotalsMedians.txt
touch readTotals.txt
touch readTotalsMedians.txt

while read mut
do

echo $mut
echo $mut >> logs/5_GetReadTotals.err


# Fill missing positions
FILE=../1_AlignData/${mut}-1_stt_pos_sort_All.txt
if [ -f "$FILE" ]; then
    echo "_All_ rep1 file exists."
else 
script=../1_AlignData/Scripts/FillMissingPositionsBedGraph.py
echo filling missing positions
python $script -p ../1_AlignData/${mut}-1_stt_pos_sort.bedGraph -m ../1_AlignData/${mut}-1_stt_neg_sort.bedGraph
echo done
fi

FILE=../1_AlignData/${mut}-2_stt_pos_sort_All.txt
if [ -f "$FILE" ]; then
    echo "_All_ rep2 file exists."
else
script=../1_AlignData/Scripts/FillMissingPositionsBedGraph.py
echo filling missing positions
python $script -p ../1_AlignData/${mut}-2_stt_pos_sort.bedGraph -m ../1_AlignData/${mut}-2_stt_neg_sort.bedGraph
echo done
fi

# Make bed files with only protein-coding regions
FILE=../1_AlignData/${mut}-1_stt_pos_sort_All_pctxpts.txt
if [ -f "$FILE" ]; then
    echo "pctxpts rep1 file exists."
else 
# rep1
bedtools intersect -a ../1_AlignData/${mut}-1_stt_pos_sort_All.txt -b ../0_Annotations/Genes/txpts_all_pos.bed > ../1_AlignData/${mut}-1_stt_pos_sort_All_pctxpts.txt
bedtools intersect -a ../1_AlignData/${mut}-1_stt_neg_sort_All.txt -b ../0_Annotations/Genes/txpts_all_neg.bed > ../1_AlignData/${mut}-1_stt_neg_sort_All_pctxpts.txt
../1_AlignData/rco1-1_stt_neg_sort_All_pctxpts.txt
fi
FILE=../1_AlignData/${mut}-2_stt_pos_sort_All_pctxpts.txt
if [ -f "$FILE" ]; then
    echo "pctxpts rep2 file exists."
else 
# rep2
bedtools intersect -a ../1_AlignData/${mut}-2_stt_pos_sort_All.txt -b ../0_Annotations/Genes/txpts_all_pos.bed > ../1_AlignData/${mut}-2_stt_pos_sort_All_pctxpts.txt
bedtools intersect -a ../1_AlignData/${mut}-2_stt_neg_sort_All.txt -b ../0_Annotations/Genes/txpts_all_neg.bed > ../1_AlignData/${mut}-2_stt_neg_sort_All_pctxpts.txt
fi

# Get total from pos strand
Pos1=`awk 'NR > 1 { sum += $4 } END { print sum }' ../1_AlignData/${mut}-1_stt_pos_sort_All_pctxpts.txt`
Pos2=`awk 'NR > 1 { sum += $4 } END { print sum }' ../1_AlignData/${mut}-2_stt_pos_sort_All_pctxpts.txt`
# Get total from neg strand
Neg1=`awk 'NR > 1 { sum += $4 } END { print sum }' ../1_AlignData/${mut}-1_stt_neg_sort_All_pctxpts.txt`
Neg2=`awk 'NR > 1 { sum += $4 } END { print sum }' ../1_AlignData/${mut}-2_stt_neg_sort_All_pctxpts.txt`
# Sum them
((sum1=$Pos1+$Neg1))
((sum2=$Pos2+$Neg2))
# Get medians
((med=($sum1+$sum2)/2))

# Write to file
echo $mut$'-1\t'$sum1 >> readTotals.txt
echo $mut$'-2\t'$sum2 >> readTotals.txt
echo $mut$'\t'$med >> readTotalsMedians.txt

done < deletionStrainsCovFiltered.txt # tmpStrains.txt # deletionStrains.txt deletionStrainsCovFiltered.txt
