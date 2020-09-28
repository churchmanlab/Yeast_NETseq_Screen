#!/bin/bash

# Script to transform un-normalized NET-seq counts into proper format for DESeq2 (gene x sample with unnormalized data)

# Written by: K. Lachance
# Date: November 27, 2018

# Use: ./makeDESeqMatrix.sh mut1 mut2 (mut3 mut4 mut5 mut6) strain

# Make wildtype matrix (used for all strains) with full NET-seq read
# Note: only need to run this once for first deletion strain analyzed
wt1_pos_file="/n/groups/churchman/kcl19/Screen/ALIGNMENT/CompletedCoverageFiles/YSC001_ES_216_cov_pos.bedGraph"
wt1_neg_file="/n/groups/churchman/kcl19/Screen/ALIGNMENT/CompletedCoverageFiles/YSC001_ES_216_cov_neg.bedGraph"
wt2_pos_file="/n/groups/churchman/kcl19/Screen/ALIGNMENT/CompletedCoverageFiles/YSC001_KH_216_cov_pos.bedGraph"
wt2_neg_file="/n/groups/churchman/kcl19/Screen/ALIGNMENT/CompletedCoverageFiles/YSC001_KH_216_cov_neg.bedGraph"
wt3_pos_file="/n/groups/churchman/kcl19/Screen/ALIGNMENT/CompletedCoverageFiles/YSC001_KT_216_cov_pos.bedGraph"
wt3_neg_file="/n/groups/churchman/kcl19/Screen/ALIGNMENT/CompletedCoverageFiles/YSC001_KT_216_cov_neg.bedGraph"
wt4_pos_file="/n/groups/churchman/kcl19/Screen/ALIGNMENT/CompletedCoverageFiles/YSC001_KT16_all_cov_pos.bedGraph"
wt4_neg_file="/n/groups/churchman/kcl19/Screen/ALIGNMENT/CompletedCoverageFiles/YSC001_KT16_all_cov_neg.bedGraph"

# Overlap individually with gene files
sortBed -i $wt1_pos_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_pos.bed - > tmp.wt1.pos.bed
sortBed -i $wt1_neg_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_neg.bed - > tmp.wt1.neg.bed
sortBed -i $wt2_pos_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_pos.bed - > tmp.wt2.pos.bed
sortBed -i $wt2_neg_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_neg.bed - > tmp.wt2.neg.bed
sortBed -i $wt3_pos_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_pos.bed - > tmp.wt3.pos.bed
sortBed -i $wt3_neg_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_neg.bed - > tmp.wt3.neg.bed
sortBed -i $wt4_pos_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_pos.bed - > tmp.wt4.pos.bed
sortBed -i $wt4_neg_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_neg.bed - > tmp.wt4.neg.bed

# Combine files by replicate
paste tmp.wt1.pos.bed tmp.wt2.pos.bed tmp.wt3.pos.bed tmp.wt4.pos.bed | awk -F "\t" 'BEGIN {OFS="\t"} {print $7, $14, $21, $28}' > tmp.wt.pos.bed
paste tmp.wt1.neg.bed tmp.wt2.neg.bed tmp.wt3.neg.bed tmp.wt4.neg.bed | awk -F "\t" 'BEGIN {OFS="\t"} {print $7, $14, $21, $28}' > tmp.wt.neg.bed

# Check number of parameters
if [ "$#" -eq 3 ]; then
    
    # If two replicates
    mut1=$1
    mut2=$2
    strain=$3
    echo "Two replicates of ${strain}: ${mut1} and ${mut2}"

    # Get input files
    mut1_pos_file="/n/groups/churchman/kcl19/Screen/ALIGNMENT/CompletedCoverageFiles/${mut1}_cov_pos.bedGraph"
    mut1_neg_file="/n/groups/churchman/kcl19/Screen/ALIGNMENT/CompletedCoverageFiles/${mut1}_cov_neg.bedGraph"
    mut2_pos_file="/n/groups/churchman/kcl19/Screen/ALIGNMENT/CompletedCoverageFiles/${mut2}_cov_pos.bedGraph"
    mut2_neg_file="/n/groups/churchman/kcl19/Screen/ALIGNMENT/CompletedCoverageFiles/${mut2}_cov_neg.bedGraph"
    
    # Overlap individually with gene files
    sortBed -i $mut1_pos_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_pos.bed - > tmp.mut1.pos.bed
    sortBed -i $mut1_neg_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_neg.bed - > tmp.mut1.neg.bed
    sortBed -i $mut2_pos_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_pos.bed - > tmp.mut2.pos.bed
    sortBed -i $mut2_neg_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_neg.bed - > tmp.mut2.neg.bed

    # Combine files by replicate
    paste tmp.mut1.pos.bed tmp.mut2.pos.bed | awk -F "\t" 'BEGIN {OFS="\t"} {print $4, $7, $14}' | paste - tmp.wt.pos.bed > tmp.pos.bed
    paste tmp.mut1.neg.bed tmp.mut2.neg.bed | awk -F "\t" 'BEGIN {OFS="\t"} {print $4, $7, $14}' | paste - tmp.wt.neg.bed > tmp.neg.bed

    # Combine files by strand
    echo -e "Gene\t${strain}-1\t${strain}-2\twt-1\twt-2\twt-3\twt-4" > header.txt

    cat header.txt tmp.pos.bed tmp.neg.bed > ${strain}_subgenic.mat.txt

    # Create conditions file
    echo -e "\tCondition\tType\n${strain}-1\ttreated\tsingle-read\n${strain}-2\ttreated\tsingle-read\nwt-1\tuntreated\tsingle-read\nwt-2\tuntreated\tsingle-read\nwt-3\tuntreated\tsingle-read\nwt-4\tuntreated\tsingle-read" > ${strain}_subgenic.cond.txt

    # Remove temporary files
    rm -fr tmp.mut*.bed tmp.pos.bed tmp.neg.bed header.txt

fi

if [ "$#" -eq 4 ]; then

    # If three replicates
    mut1=$1
    mut2=$2
    mut3=$3
    strain=$4
    echo "Three replicates of ${strain}: ${mut1}, ${mut2}, and ${mut3}"

    # Get input files
    mut1_pos_file="../ALIGNMENT/CompletedCoverageFiles/${mut1}_cov_pos.bedGraph"
    mut1_neg_file="../ALIGNMENT/CompletedCoverageFiles/${mut1}_cov_neg.bedGraph"
    mut2_pos_file="../ALIGNMENT/CompletedCoverageFiles/${mut2}_cov_pos.bedGraph"
    mut2_neg_file="../ALIGNMENT/CompletedCoverageFiles/${mut2}_cov_neg.bedGraph"
    mut3_pos_file="../ALIGNMENT/CompletedCoverageFiles/${mut3}_cov_pos.bedGraph"
    mut3_neg_file="../ALIGNMENT/CompletedCoverageFiles/${mut3}_cov_neg.bedGraph"

    # Overlap individually with gene files
    sortBed -i $mut1_pos_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_pos.bed - > tmp.mut1.pos.bed
    sortBed -i $mut1_neg_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_neg.bed - > tmp.mut1.neg.bed
    sortBed -i $mut2_pos_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_pos.bed - > tmp.mut2.pos.bed
    sortBed -i $mut2_neg_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_neg.bed - > tmp.mut2.neg.bed
    sortBed -i $mut3_pos_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_pos.bed - > tmp.mut3.pos.bed
    sortBed -i $mut3_neg_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_neg.bed - > tmp.mut3.neg.bed

    # Combine files by replicate
    paste tmp.mut1.pos.bed tmp.mut2.pos.bed tmp.mut3.pos.bed | awk -F "\t" 'BEGIN {OFS="\t"} {print $4, $7, $14, $21}' | paste - tmp.wt.pos.bed > tmp.pos.bed
    paste tmp.mut1.neg.bed tmp.mut2.neg.bed tmp.mut3.neg.bed | awk -F "\t" 'BEGIN {OFS="\t"} {print $4, $7, $14, $21}' | paste - tmp.wt.neg.bed > tmp.neg.bed

    # Combine files by strand
    echo -e "Gene\t${strain}-1\t${strain}-2\t${strain}-3\twt-1\twt-2\twt-3\twt-4" > header.txt

    cat header.txt tmp.pos.bed tmp.neg.bed > ${strain}_subgenic.mat.txt

    # Create conditions file
    echo -e "\tCondition\tType\n${strain}-1\ttreated\tsingle-read\n${strain}-2\ttreated\tsingle-read\n${strain}-3\ttreated\tsingle-read\nwt-1\tuntreated\tsingle-read\nwt-2\tuntreated\tsingle-read\nwt-3\tuntreated\tsingle-read\nwt-4\tuntreated\tsingle-read" > ${strain}_subgenic.cond.txt

    # Remove temporary files
    rm -fr tmp.mut*.bed tmp.pos.bed tmp.neg.bed header.txt

fi

if [ "$#" -eq 5 ]; then

    # If four replicates
    mut1=$1
    mut2=$2
    mut3=$3
    mut4=$4
    strain=$5
    echo "Four replicates of ${strain}: ${mut1}, ${mut2}, ${mut3}, and ${mut4}"

    # Get input files
    mut1_pos_file="../ALIGNMENT/CompletedCoverageFiles/${mut1}_cov_pos.bedGraph"
    mut1_neg_file="../ALIGNMENT/CompletedCoverageFiles/${mut1}_cov_neg.bedGraph"
    mut2_pos_file="../ALIGNMENT/CompletedCoverageFiles/${mut2}_cov_pos.bedGraph"
    mut2_neg_file="../ALIGNMENT/CompletedCoverageFiles/${mut2}_cov_neg.bedGraph"
    mut3_pos_file="../ALIGNMENT/CompletedCoverageFiles/${mut3}_cov_pos.bedGraph"
    mut3_neg_file="../ALIGNMENT/CompletedCoverageFiles/${mut3}_cov_neg.bedGraph"
    mut4_pos_file="../ALIGNMENT/CompletedCoverageFiles/${mut4}_cov_pos.bedGraph"
    mut4_neg_file="../ALIGNMENT/CompletedCoverageFiles/${mut4}_cov_neg.bedGraph"

    # Overlap individually with gene files
    sortBed -i $mut1_pos_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_pos.bed - > tmp.mut1.pos.bed
    sortBed -i $mut1_neg_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_neg.bed - > tmp.mut1.neg.bed
    sortBed -i $mut2_pos_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_pos.bed - > tmp.mut2.pos.bed
    sortBed -i $mut2_neg_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_neg.bed - > tmp.mut2.neg.bed
    sortBed -i $mut3_pos_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_pos.bed - > tmp.mut3.pos.bed
    sortBed -i $mut3_neg_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_neg.bed - > tmp.mut3.neg.bed
    sortBed -i $mut4_pos_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_pos.bed - > tmp.mut4.pos.bed
    sortBed -i $mut4_neg_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_neg.bed - > tmp.mut4.neg.bed

    # Combine files by replicate
    paste tmp.mut1.pos.bed tmp.mut2.pos.bed tmp.mut3.pos.bed tmp.mut4.pos.bed | awk -F "\t" 'BEGIN {OFS="\t"} {print $4, $7, $14, $21, $28}' | paste - tmp.wt.pos.bed > tmp.pos.bed
    paste tmp.mut1.neg.bed tmp.mut2.neg.bed tmp.mut3.neg.bed tmp.mut4.neg.bed | awk -F "\t" 'BEGIN {OFS="\t"} {print $4, $7, $14, $21, $28}' | paste - tmp.wt.neg.bed > tmp.neg.bed

    # Combine files by strand
    echo -e "Gene\t${strain}-1\t${strain}-2\t${strain}-3\t${strain}-4\twt-1\twt-2\twt-3\twt-4" > header.txt

    cat header.txt tmp.pos.bed tmp.neg.bed > ${strain}_subgenic.mat.txt

    # Create conditions file
    echo -e "\tCondition\tType\n${strain}-1\ttreated\tsingle-read\n${strain}-2\ttreated\tsingle-read\n${strain}-3\ttreated\tsingle-read\n${strain}-4\ttreated\tsingle-read\nwt-1\tuntreated\tsingle-read\nwt-2\tuntreated\tsingle-read\nwt-3\tuntreated\tsingle-read\nwt-4\tuntreated\tsingle-read" > ${strain}_subgenic.cond.txt

    # Remove temporary files
    rm -fr tmp.mut*.bed tmp.pos.bed tmp.neg.bed header.txt

fi

if [ "$#" -eq 6 ]; then

    # If five replicates
    mut1=$1
    mut2=$2
    mut3=$3
    mut4=$4
    mut5=$5
    strain=$6
    echo "Five replicates of ${strain}: ${mut1}, ${mut2}, ${mut3}, ${mut4}, and ${mut5}"

    # Get input files
    mut1_pos_file="../ALIGNMENT/CompletedCoverageFiles/${mut1}_cov_pos.bedGraph"
    mut1_neg_file="../ALIGNMENT/CompletedCoverageFiles/${mut1}_cov_neg.bedGraph"
    mut2_pos_file="../ALIGNMENT/CompletedCoverageFiles/${mut2}_cov_pos.bedGraph"
    mut2_neg_file="../ALIGNMENT/CompletedCoverageFiles/${mut2}_cov_neg.bedGraph"
    mut3_pos_file="../ALIGNMENT/CompletedCoverageFiles/${mut3}_cov_pos.bedGraph"
    mut3_neg_file="../ALIGNMENT/CompletedCoverageFiles/${mut3}_cov_neg.bedGraph"
    mut4_pos_file="../ALIGNMENT/CompletedCoverageFiles/${mut4}_cov_pos.bedGraph"
    mut4_neg_file="../ALIGNMENT/CompletedCoverageFiles/${mut4}_cov_neg.bedGraph"
    mut5_pos_file="../ALIGNMENT/CompletedCoverageFiles/${mut5}_cov_pos.bedGraph"
    mut5_neg_file="../ALIGNMENT/CompletedCoverageFiles/${mut5}_cov_neg.bedGraph"

    # Overlap individually with gene files
    sortBed -i $mut1_pos_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_pos.bed - > tmp.mut1.pos.bed
    sortBed -i $mut1_neg_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_neg.bed - > tmp.mut1.neg.bed
    sortBed -i $mut2_pos_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_pos.bed - > tmp.mut2.pos.bed
    sortBed -i $mut2_neg_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_neg.bed - > tmp.mut2.neg.bed
    sortBed -i $mut3_pos_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_pos.bed - > tmp.mut3.pos.bed
    sortBed -i $mut3_neg_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_neg.bed - > tmp.mut3.neg.bed
    sortBed -i $mut4_pos_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_pos.bed - > tmp.mut4.pos.bed
    sortBed -i $mut4_neg_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_neg.bed - > tmp.mut4.neg.bed
    sortBed -i $mut5_pos_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_pos.bed - > tmp.mut5.pos.bed
    sortBed -i $mut5_neg_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_neg.bed - > tmp.mut5.neg.bed

    # Combine files by replicate
    paste tmp.mut1.pos.bed tmp.mut2.pos.bed tmp.mut3.pos.bed tmp.mut4.pos.bed tmp.mut5.pos.bed| awk -F "\t" 'BEGIN {OFS="\t"} {print $4, $7, $14, $21, $28, $35}' | paste - tmp.wt.pos.bed > tmp.pos.bed
    paste tmp.mut1.neg.bed tmp.mut2.neg.bed tmp.mut3.neg.bed tmp.mut4.neg.bed tmp.mut5.neg.bed | awk -F "\t" 'BEGIN {OFS="\t"} {print $4, $7, $14, $21, $28, $35}'| paste - tmp.wt.neg.bed > tmp.neg.bed

    # Combine files by strand
    echo -e "Gene\t${strain}-1\t${strain}-2\t${strain}-3\t${strain}-4\t${strain}-5\twt-1\twt-2\twt-3\twt-4" > header.txt

    cat header.txt tmp.pos.bed tmp.neg.bed > ${strain}_subgenic.mat.txt

    # Create conditions file
    echo -e "\tCondition\tType\n${strain}-1\ttreated\tsingle-read\n${strain}-2\ttreated\tsingle-read\n${strain}-3\ttreated\tsingle-read\n${strain}-4\ttreated\tsingle-read\n${strain}-5\ttreated\tsingle-read\nwt-1\tuntreated\tsingle-read\nwt-2\tuntreated\tsingle-read\nwt-3\tuntreated\tsingle-read\nwt-4\tuntreated\tsingle-read" > ${strain}_subgenic.cond.txt

    # Remove temporary files
    rm -fr tmp.mut*.bed tmp.pos.bed tmp.neg.bed header.txt

fi

if [ "$#" -eq 7 ]; then

    # If six replicates
    mut1=$1
    mut2=$2
    mut3=$3
    mut4=$4
    mut5=$5
    mut6=$6
    strain=$7
    echo "Six replicates of ${strain}: ${mut1}, ${mut2}, ${mut3}, ${mut4}, ${mut5}, and ${mut6}"

    # Get input files
    mut1_pos_file="../ALIGNMENT/CompletedCoverageFiles/${mut1}_cov_pos.bedGraph"
    mut1_neg_file="../ALIGNMENT/CompletedCoverageFiles/${mut1}_cov_neg.bedGraph"
    mut2_pos_file="../ALIGNMENT/CompletedCoverageFiles/${mut2}_cov_pos.bedGraph"
    mut2_neg_file="../ALIGNMENT/CompletedCoverageFiles/${mut2}_cov_neg.bedGraph"
    mut3_pos_file="../ALIGNMENT/CompletedCoverageFiles/${mut3}_cov_pos.bedGraph"
    mut3_neg_file="../ALIGNMENT/CompletedCoverageFiles/${mut3}_cov_neg.bedGraph"
    mut4_pos_file="../ALIGNMENT/CompletedCoverageFiles/${mut4}_cov_pos.bedGraph"
    mut4_neg_file="../ALIGNMENT/CompletedCoverageFiles/${mut4}_cov_neg.bedGraph"
    mut5_pos_file="../ALIGNMENT/CompletedCoverageFiles/${mut5}_cov_pos.bedGraph"
    mut5_neg_file="../ALIGNMENT/CompletedCoverageFiles/${mut5}_cov_neg.bedGraph"
    mut6_pos_file="../ALIGNMENT/CompletedCoverageFiles/${mut6}_cov_pos.bedGraph"
    mut6_neg_file="../ALIGNMENT/CompletedCoverageFiles/${mut6}_cov_neg.bedGraph"

    # Overlap individually with gene files
    sortBed -i $mut1_pos_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_pos.bed - > tmp.mut1.pos.bed
    sortBed -i $mut1_neg_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_neg.bed - > tmp.mut1.neg.bed
    sortBed -i $mut2_pos_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_pos.bed - > tmp.mut2.pos.bed
    sortBed -i $mut2_neg_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_neg.bed - > tmp.mut2.neg.bed
    sortBed -i $mut3_pos_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_pos.bed - > tmp.mut3.pos.bed
    sortBed -i $mut3_neg_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_neg.bed - > tmp.mut3.neg.bed
    sortBed -i $mut4_pos_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_pos.bed - > tmp.mut4.pos.bed
    sortBed -i $mut4_neg_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_neg.bed - > tmp.mut4.neg.bed
    sortBed -i $mut5_pos_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_pos.bed - > tmp.mut5.pos.bed
    sortBed -i $mut5_neg_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_neg.bed - > tmp.mut5.neg.bed
    sortBed -i $mut6_pos_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_pos.bed - > tmp.mut6.pos.bed
    sortBed -i $mut6_neg_file | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, ".", $4}' - | bedmap --echo --delim "\t" --sum ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED+150-100_neg.bed - > tmp.mut6.neg.bed

    # Combine files by replicate
    paste tmp.mut1.pos.bed tmp.mut2.pos.bed tmp.mut3.pos.bed tmp.mut4.pos.bed tmp.mut5.pos.bed tmp.mut6.pos.bed | awk -F "\t" 'BEGIN {OFS="\t"} {print $4, $7, $14, $21, $28, $35, $42}'| paste - tmp.wt.pos.bed > tmp.pos.bed
    paste tmp.mut1.neg.bed tmp.mut2.neg.bed tmp.mut3.neg.bed tmp.mut4.neg.bed tmp.mut5.neg.bed tmp.mut6.neg.bed | awk -F "\t" 'BEGIN {OFS="\t"} {print $4, $7, $14, $21, $28, $35, $42}'| paste - tmp.wt.neg.bed > tmp.neg.bed

    # Combine files by strand
    echo -e "Gene\t${strain}-1\t${strain}-2\t${strain}-3\t${strain}-4\t${strain}-5\t${strain}-6\twt-1\twt-2\twt-3\twt-4" > header.txt

    cat header.txt tmp.pos.bed tmp.neg.bed > ${strain}_subgenic.mat.txt

    # Create conditions file
    echo -e "\tCondition\tType\n${strain}-1\ttreated\tsingle-read\n${strain}-2\ttreated\tsingle-read\n${strain}-3\ttreated\tsingle-read\n${strain}-4\ttreated\tsingle-read\n${strain}-5\ttreated\tsingle-read\n${strain}-6\ttreated\tsingle-read\nwt-1\tuntreated\tsingle-read\nwt-2\tuntreated\tsingle-read\nwt-3\tuntreated\tsingle-read\nwt-4\tuntreated\tsingle-read" > ${strain}_subgenic.cond.txt

    # Remove temporary files
    rm -fr tmp.mut*.bed tmp.pos.bed tmp.neg.bed header.txt

fi

# Remove temporary files created from wildtype analysis (uncomment when running this script for last time)
rm -fr tmp.wt*.bed
