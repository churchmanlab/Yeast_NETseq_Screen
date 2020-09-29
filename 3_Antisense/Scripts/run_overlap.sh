#!/bin/bash

while read mut
do

    echo $mut

    # Note swap of strands to look at antisense 
    ./Scripts/overlap_data_whole.R /n/groups/churchman/kcl19/Screen/DATA/CombinedReplicates/${mut}.pos.bed ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED-250+4000_neg.bed neg ${mut}
    ./Scripts/overlap_data_whole.R /n/groups/churchman/kcl19/Screen/DATA/CombinedReplicates/${mut}.neg.bed ../0_Annotations/Genes/genes_polII.noOver.coding.TRIMMED-250+4000_pos.bed pos ${mut}

    # Clean up temporary files
    rm -fr *_tmpCoord_*

    # If using SLURM scheduler
    sbatch -p short -t 0-1:00 --job-name $mut --mem 50G -o ${mut}.out -e ${mut}.err --wrap ""./Scripts/plotWholeGeneHeatmap.R" "${mut}"" 
    
    # If not using SLURM scheduler
    ./Scripts/plotWholeGeneHeatmap.R ${mut}

done < deletionStrains.txt
