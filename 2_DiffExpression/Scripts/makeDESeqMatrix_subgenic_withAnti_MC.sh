#!/bin/bash

#SBATCH -c 1
#SBATCH -t 0-01:00
#SBATCH -p short
#SBATCH --mem=5G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=youremail@gmail.com

### Script to transform un-normalized NET-seq counts into proper format for DESeq2 (gene x sample with unnormalized data)

# Written by: K. Lachance
# Date: November 27, 2018
# Updated 12/2021 by M. Couvillion to use gene file from Xu et al that also includes antisense transcripts that I added in

# Use: ./makeDESeqMatrix_subgenic_withAnti_MC.sh mut1 mut2 (mut3 mut4 mut5 mut6) strain

# Run first with one strain then comment out wt lines and repeat with rest
# muts="bre1"
# muts="bur2 bye1 cac1 cac2 cac3 cbc1 ccr4 cdc39 cdc73 chd1 ctr9 dhh1 dst1 eaf1 elf1 gcn5 hda1 hpc2 htz1 ino80 isw1 isw2 leo1 nap1 npl3 paf1 rad6 rco1 rpb4 rph1 rsc30 rtf1 rtt103 set2 set3 spt4 swr1 ubp8 vps15 vps34"
# for mut in $muts
# do
# sbatch -o logs/makeMatrix_subgenic_${mut}.log -e logs/makeMatrix_subgenic_${mut}.err ../Scripts/makeDESeqMatrix_subgenic_withAnti_MC.sh ${mut}-1 ${mut}-2 $mut
# done
# 

module load bedops/2.4.30


pos_gene_file_in="/n/groups/churchman/mc348/NETseqScreen/Yeast_NETseq_Screen-master/0_Annotations/Genes/txpts_all_and_antisense_pos_sort_mergeByGeneTRIMMEDp150m100.bed"
neg_gene_file_in="/n/groups/churchman/mc348/NETseqScreen/Yeast_NETseq_Screen-master/0_Annotations/Genes/txpts_all_and_antisense_neg_sort_mergeByGeneTRIMMEDp150m100.bed"

# Have to resort these to use bedmap
sort-bed $pos_gene_file_in > ${pos_gene_file_in/.bed/_bedopsort.bed}
sort-bed $neg_gene_file_in > ${neg_gene_file_in/.bed/_bedopsort.bed}

pos_gene_file=${pos_gene_file_in/.bed/_bedopsort.bed}
neg_gene_file=${neg_gene_file_in/.bed/_bedopsort.bed}

# Make wildtype matrix (used for all strains) with full NET-seq read
# Note: only need to run this once for first deletion strain analyzed
wt1_pos_file="/n/groups/churchman/mc348/NETseqScreen/Yeast_NETseq_Screen-master/1_AlignData/wt-1_stt_pos_sort.bedGraph" # "../1_AlignData/wt-1_cov_pos.bedGraph"
wt1_neg_file="/n/groups/churchman/mc348/NETseqScreen/Yeast_NETseq_Screen-master/1_AlignData/wt-1_stt_neg_sort.bedGraph"
wt2_pos_file="/n/groups/churchman/mc348/NETseqScreen/Yeast_NETseq_Screen-master/1_AlignData/wt-2_stt_pos_sort.bedGraph"
wt2_neg_file="/n/groups/churchman/mc348/NETseqScreen/Yeast_NETseq_Screen-master/1_AlignData/wt-2_stt_neg_sort.bedGraph"
wt3_pos_file="/n/groups/churchman/mc348/NETseqScreen/Yeast_NETseq_Screen-master/1_AlignData/wt-3_stt_pos_sort.bedGraph"
wt3_neg_file="/n/groups/churchman/mc348/NETseqScreen/Yeast_NETseq_Screen-master/1_AlignData/wt-3_stt_neg_sort.bedGraph"
wt4_pos_file="/n/groups/churchman/mc348/NETseqScreen/Yeast_NETseq_Screen-master/1_AlignData/wt-4_stt_pos_sort.bedGraph"
wt4_neg_file="/n/groups/churchman/mc348/NETseqScreen/Yeast_NETseq_Screen-master/1_AlignData/wt-4_stt_neg_sort.bedGraph"


# Overlap individually with gene files 
sed '1d' $wt1_pos_file | sort-bed - | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum $pos_gene_file - > tmp.sub_wt1.pos.bed
sed '1d' $wt1_neg_file | sort-bed - | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum $neg_gene_file - > tmp.sub_wt1.neg.bed
sed '1d' $wt2_pos_file | sort-bed - | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum $pos_gene_file - > tmp.sub_wt2.pos.bed
sed '1d' $wt2_neg_file | sort-bed - | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum $neg_gene_file - > tmp.sub_wt2.neg.bed
sed '1d' $wt3_pos_file | sort-bed - | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum $pos_gene_file - > tmp.sub_wt3.pos.bed
sed '1d' $wt3_neg_file | sort-bed - | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum $neg_gene_file - > tmp.sub_wt3.neg.bed
sed '1d' $wt4_pos_file | sort-bed - | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum $pos_gene_file - > tmp.sub_wt4.pos.bed
sed '1d' $wt4_neg_file | sort-bed - | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum $neg_gene_file - > tmp.sub_wt4.neg.bed

# Combine files by replicate
paste tmp.sub_wt1.pos.bed tmp.sub_wt2.pos.bed tmp.sub_wt3.pos.bed tmp.sub_wt4.pos.bed | awk -F "\t" 'BEGIN {OFS="\t"} {print $7, $14, $21, $28}' > wt_sub.pos.counts.bed
paste tmp.sub_wt1.neg.bed tmp.sub_wt2.neg.bed tmp.sub_wt3.neg.bed tmp.sub_wt4.neg.bed | awk -F "\t" 'BEGIN {OFS="\t"} {print $7, $14, $21, $28}' > wt_sub.neg.counts.bed

# Check number of parameters
if [ "$#" -eq 3 ]; then
    
    # If two replicates
    mut1=$1
    mut2=$2
    strain=$3
    echo "Two replicates of ${strain}: ${mut1} and ${mut2}"

    # Get input files
    mut1_pos_file="/n/groups/churchman/mc348/NETseqScreen/Yeast_NETseq_Screen-master/1_AlignData/${mut1}_stt_pos_sort.bedGraph"
    mut1_neg_file="/n/groups/churchman/mc348/NETseqScreen/Yeast_NETseq_Screen-master/1_AlignData/${mut1}_stt_neg_sort.bedGraph"
    mut2_pos_file="/n/groups/churchman/mc348/NETseqScreen/Yeast_NETseq_Screen-master/1_AlignData/${mut2}_stt_pos_sort.bedGraph"
    mut2_neg_file="/n/groups/churchman/mc348/NETseqScreen/Yeast_NETseq_Screen-master/1_AlignData/${mut2}_stt_neg_sort.bedGraph"
    
    # Overlap individually with gene files
    sed '1d' $mut1_pos_file | sort-bed - | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum $pos_gene_file - > tmp.${mut1}.pos.bed
    sed '1d' $mut1_neg_file | sort-bed - | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum $neg_gene_file - > tmp.${mut1}.neg.bed
    sed '1d' $mut2_pos_file | sort-bed - | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum $pos_gene_file - > tmp.${mut2}.pos.bed
    sed '1d' $mut2_neg_file | sort-bed - | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum $neg_gene_file - > tmp.${mut2}.neg.bed

    # Combine files by replicate
    paste tmp.${mut1}.pos.bed tmp.${mut2}.pos.bed | awk -F "\t" 'BEGIN {OFS="\t"} {print $4, $7, $14}' | paste - wt_sub.pos.counts.bed > tmp.${strain}.pos.bed
    paste tmp.${mut1}.neg.bed tmp.${mut2}.neg.bed | awk -F "\t" 'BEGIN {OFS="\t"} {print $4, $7, $14}' | paste - wt_sub.neg.counts.bed > tmp.${strain}.neg.bed

    # Combine files by strand
    echo -e "Gene\t${strain}-1\t${strain}-2\twt-1\twt-2\twt-3\twt-4" > ${strain}_sub_header.txt

    cat ${strain}_sub_header.txt tmp.${strain}.pos.bed tmp.${strain}.neg.bed > ${strain}_subgenic.mat.txt
    

    # Create conditions file
    echo -e "\tCondition\tType\n${strain}-1\ttreated\tsingle-read\n${strain}-2\ttreated\tsingle-read\nwt-1\tuntreated\tsingle-read\nwt-2\tuntreated\tsingle-read\nwt-3\tuntreated\tsingle-read\nwt-4\tuntreated\tsingle-read" > ${strain}_subgenic.cond.txt

    # Remove temporary files
    rm -fr tmp.${strain}*.bed ${strain}_sub_header.txt

fi


# Remove temporary files created from wildtype analysis (uncomment when running this script for last time)
# rm -fr tmp.wt*.bed
