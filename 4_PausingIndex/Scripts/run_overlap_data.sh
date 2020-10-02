#!/bin/bash

while read mut
do

    echo $mut

#    ./Scripts/overlap_data.R /n/groups/churchman/kcl19/Screen/DATA/CombinedReplicates/${mut}.pos.bed ../0_Annotations/Genes/genesTrimmedTSS100-600.pos.bed pos ${mut}.TSS
#    ./Scripts/overlap_data.R /n/groups/churchman/kcl19/Screen/DATA/CombinedReplicates/${mut}.neg.bed  ../0_Annotations/Genes/genesTrimmedTSS100-600.neg.bed neg ${mut}.TSS

#    ./Scripts/overlap_data.R /n/groups/churchman/kcl19/Screen/DATA/CombinedReplicates/${mut}.pos.bed ../0_Annotations/Genes/genesTrimmedpA500-200.pos.bed pos ${mut}.pA
#    ./Scripts/overlap_data.R /n/groups/churchman/kcl19/Screen/DATA/CombinedReplicates/${mut}.neg.bed ../0_Annotations/Genes/genesTrimmedpA500-200.neg.bed neg ${mut}.pA

#    ./Scripts/overlap_data.R /n/groups/churchman/kcl19/Screen/DATA/CombinedReplicates/${mut}.pos.bed ../0_Annotations/Genes/genesTrimmedpA500-200.neg.bed neg ${mut}.anti
#    ./Scripts/overlap_data.R /n/groups/churchman/kcl19/Screen/DATA/CombinedReplicates/${mut}.neg.bed ../0_Annotations/Genes/genesTrimmedpA500-200.pos.bed pos ${mut}.anti

     ./Scripts/overlap_data.R /n/groups/churchman/kcl19/Screen/DATA/CombinedReplicates/${mut}.pos.bed ../0_Annotations/Genes/genesTrimmed5pSS25.pos.bed pos ${mut}.5pSS
     ./Scripts/overlap_data.R /n/groups/churchman/kcl19/Screen/DATA/CombinedReplicates/${mut}.neg.bed ../0_Annotations/Genes/genesTrimmed5pSS25.neg.bed neg ${mut}.5pSS
     
     ./Scripts/overlap_data.R /n/groups/churchman/kcl19/Screen/DATA/CombinedReplicates/${mut}.pos.bed ../0_Annotations/Genes/genesTrimmed3pSS25.pos.bed pos ${mut}.3pSS
     ./Scripts/overlap_data.R /n/groups/churchman/kcl19/Screen/DATA/CombinedReplicates/${mut}.neg.bed ../0_Annotations/Genes/genesTrimmed3pSS25.neg.bed neg ${mut}.3pSS



done < tmp
#done < deletionStrains.txt

# Remove temporary files
rm -fr *_tmpCoord_*
