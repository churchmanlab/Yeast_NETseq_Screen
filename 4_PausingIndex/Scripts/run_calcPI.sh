#!/bin/bash

while read mut
do

    echo $mut
    
    # TSSs
    echo "     Working on TSSs..."

    intersectBed -a ../0_Annotations/Genes/genesTrimmedTSS.pos.bed -b /n/groups/churchman/kcl19/Screen/DATA/CombinedReplicates/${mut}.noNorm.pos.bed -wa -wb > ${mut}.tmp.pos.bed
    intersectBed -a ../0_Annotations/Genes/genesTrimmedTSS.neg.bed -b /n/groups/churchman/kcl19/Screen/DATA/CombinedReplicates/${mut}.noNorm.neg.bed -wa -wb > ${mut}.tmp.neg.bed

    intersectBed -a ../0_Annotations/Genes/genesTrimmedTSS_NOT.pos.bed -b /n/groups/churchman/kcl19/Screen/DATA/CombinedReplicates/${mut}.noNorm.pos.bed -wa -wb > ${mut}.tmpNOT.pos.bed
    intersectBed -a ../0_Annotations/Genes/genesTrimmedTSS_NOT.neg.bed -b /n/groups/churchman/kcl19/Screen/DATA/CombinedReplicates/${mut}.noNorm.neg.bed -wa -wb > ${mut}.tmpNOT.neg.bed

    ./Scripts/calcPI.R ${mut} TSS
    rm -fr ${mut}.tmp*

    # pA sites
    echo "     Working on pA sites..."

    intersectBed -a ../0_Annotations/Genes/genesTrimmedpA.pos.bed -b /n/groups/churchman/kcl19/Screen/DATA/CombinedReplicates/${mut}.noNorm.pos.bed -wa -wb > ${mut}.tmp.pos.bed
    intersectBed -a ../0_Annotations/Genes/genesTrimmedpA.neg.bed -b /n/groups/churchman/kcl19/Screen/DATA/CombinedReplicates/${mut}.noNorm.neg.bed -wa -wb > ${mut}.tmp.neg.bed

    intersectBed -a ../0_Annotations/Genes/genesTrimmedpA_NOT.pos.bed -b /n/groups/churchman/kcl19/Screen/DATA/CombinedReplicates/${mut}.noNorm.pos.bed -wa -wb > ${mut}.tmpNOT.pos.bed
    intersectBed -a ../0_Annotations/Genes/genesTrimmedpA_NOT.neg.bed -b /n/groups/churchman/kcl19/Screen/DATA/CombinedReplicates/${mut}.noNorm.neg.bed -wa -wb > ${mut}.tmpNOT.neg.bed

    ./Scripts/calcPI.R ${mut} pA
    rm -fr ${mut}.tmp*

    # 5' splice sites
    echo "     Working on 5' splice sites..."

    intersectBed -a ../0_Annotations/Genes/genesTrimmed5pSS.pos.bed -b /n/groups/churchman/kcl19/Screen/DATA/CombinedReplicates/${mut}.noNorm.pos.bed -wa -wb > ${mut}.tmp.pos.bed
    intersectBed -a ../0_Annotations/Genes/genesTrimmed5pSS.neg.bed -b /n/groups/churchman/kcl19/Screen/DATA/CombinedReplicates/${mut}.noNorm.neg.bed -wa -wb > ${mut}.tmp.neg.bed

    intersectBed -a ../0_Annotations/Genes/genesTrimmedSS_NOT.pos.bed -b /n/groups/churchman/kcl19/Screen/DATA/CombinedReplicates/${mut}.noNorm.pos.bed -wa -wb > ${mut}.tmpNOT.pos.bed
    intersectBed -a ../0_Annotations/Genes/genesTrimmedSS_NOT.neg.bed -b /n/groups/churchman/kcl19/Screen/DATA/CombinedReplicates/${mut}.noNorm.neg.bed -wa -wb > ${mut}.tmpNOT.neg.bed

    ./Scripts/calcPI.R ${mut} 5pSS
    rm -fr ${mut}.tmp*

    # 3' splice sites
    echo "     Working on 3' splice sites..."

    intersectBed -a ../0_Annotations/Genes/genesTrimmed3pSS.pos.bed -b /n/groups/churchman/kcl19/Screen/DATA/CombinedReplicates/${mut}.noNorm.pos.bed -wa -wb > ${mut}.tmp.pos.bed
    intersectBed -a ../0_Annotations/Genes/genesTrimmed3pSS.neg.bed -b /n/groups/churchman/kcl19/Screen/DATA/CombinedReplicates/${mut}.noNorm.neg.bed -wa -wb > ${mut}.tmp.neg.bed

    intersectBed -a ../0_Annotations/Genes/genesTrimmedSS_NOT.pos.bed -b /n/groups/churchman/kcl19/Screen/DATA/CombinedReplicates/${mut}.noNorm.pos.bed -wa -wb > ${mut}.tmpNOT.pos.bed
    intersectBed -a ../0_Annotations/Genes/genesTrimmedSS_NOT.neg.bed -b /n/groups/churchman/kcl19/Screen/DATA/CombinedReplicates/${mut}.noNorm.neg.bed -wa -wb > ${mut}.tmpNOT.neg.bed

    ./Scripts/calcPI.R ${mut} 3pSS
    rm -fr ${mut}.tmp*

    # Antisense
    echo "     Working on 3' antisense..."

    cat /n/groups/churchman/kcl19/Screen/DATA/CombinedReplicates/${mut}.pos.bed | awk -F "\t" '{print $1, $2, $3, ".", $4}' | sed 's/ /\t/g' | bedmap --sum ../0_Annotations/Genes/genesTrimmedAntisense.pos.bed - > tmp.5p_pos.txt
    cat /n/groups/churchman/kcl19/Screen/DATA/CombinedReplicates/${mut}.pos.bed | awk -F "\t" '{print $1, $2, $3, ".", $4}' | sed 's/ /\t/g' | bedmap --sum ../0_Annotations/Genes/genesTrimmedAntisense_NOT.pos.bed - > tmp.gb_pos.txt
    paste ../0_Annotations/Genes/genesTrimmedAntisense.pos.bed tmp.5p_pos.txt tmp.gb_pos.txt | awk -F "\t" '{if ($7 + 0 != 0) print $0, $8/$7; else print $0, "NA"}' |sed 's/ /\t/g' | awk -v mut=$mut -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, $2, $3, $4, $5, $6, "0", "0", $7, $8, $9, "Anti", mut}' > ${mut}_Antisense.PI.pos.bed

    cat /n/groups/churchman/kcl19/Screen/DATA/CombinedReplicates/${mut}.neg.bed | awk -F "\t" '{print $1, $2, $3, ".", $4}' | sed 's/ /\t/g' | bedmap --sum ../0_Annotations/Genes/genesTrimmedAntisense.neg.bed - > tmp.5p_neg.txt
    cat /n/groups/churchman/kcl19/Screen/DATA/CombinedReplicates/${mut}.neg.bed | awk -F "\t" '{print $1, $2, $3, ".", $4}' | sed 's/ /\t/g' | bedmap --sum ../0_Annotations/Genes/genesTrimmedAntisense_NOT.neg.bed - > tmp.gb_neg.txt
    paste ../0_Annotations/Genes/genesTrimmedAntisense.neg.bed tmp.5p_neg.txt tmp.gb_neg.txt | awk -F "\t" '{if ($7 + 0 != 0) print $0, $8/$7; else print $0, "NA"}' |sed 's/ /\t/g' | awk -v mut=$mut -F "\t" 'BEGIN {OFS="\t"} {print $1, $2, $3, $2, $3, $4, $5, $6, "0", "0", $7, $8, $9, "Anti", mut}' > ${mut}_Antisense.PI.neg.bed

    cat ${mut}_Antisense.PI.pos.bed ${mut}_Antisense.PI.neg.bed > ${mut}_Antisense.PI.bed

    rm -fr tmp.*.txt ${mut}_Antisense.PI.pos.bed ${mut}_Antisense.PI.neg.bed
    
#done < deletionStrains.txt
done < $1
