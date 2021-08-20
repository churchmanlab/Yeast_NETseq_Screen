#!/bin/bash

#SBATCH -c 1
#SBATCH -t 0-01:00
#SBATCH -p short
#SBATCH --mem=10G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<youremail@address.com>

### Find pauses above IDR threshold for each mut

# Use: 
# echo -e 'yfg1\ngoi1\nwt' > deletionStrains.txt
# geneCov=2
# sbatch -o logs/5_IDR_${geneCov}.log -e logs/5_IDR_${geneCov}.err ./Scripts/IDR.sh $geneCov

w=100 # half window size
t=$1 # gene coverage threshold

while read mut
do

    echo $mut

    # Get list of pauses with scores
    ./Scripts/findOverlap_pauseScore_MC.R $mut $w $t

    # Calculate and plot IDR
    ./Scripts/plot_pauseScore_MC.R $mut $t


done < deletionStrains.txt



###### Do this after all jobs complete
# To make bedGraph-like files from SampleName_PASS_pauseScore.txt
# Libs='wt'
# for lib in $Libs
# do
# awk -F'\t' 'BEGIN{OFS="\t"; print "track type=bedGraph";} { if ($6 == "+") print $1"\t"$2"\t"$3"\t10\t.\t"$6}' ${lib}_PASS_pauseScore_cov2.txt > ${lib}_PASS_pause.pos.w100.cov2.bedGraph
# awk -F'\t' 'BEGIN{OFS="\t"; print "track type=bedGraph";} { if ($6 == "-") print $1"\t"$2"\t"$3"\t10\t.\t"$6}' ${lib}_PASS_pauseScore_cov2.txt > ${lib}_PASS_pause.neg.w100.cov2.bedGraph
# cat ${lib}_PASS_pause.pos.w100.cov2.bedGraph ${lib}_PASS_pause.neg.w100.cov2.bedGraph > ${lib}_PASS_pause.w100.cov2.bedGraph
# done