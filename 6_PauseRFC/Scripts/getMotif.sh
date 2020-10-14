#!/bin/bash


# Must load all in this order to work
module load gcc/6.2.0
module load openmpi/3.1.0
module load meme/5.0.3

while read mut
do

    echo $mut
    mkdir -p ${mut}

    # Get motifs
    # NOTE: must not be on interactive mode for this to work
    meme ../5_PauseLoci/${mut}.IDRrep.pauseWindows.fasta -dna -oc . -nostatus -time 18000 -mod zoops -nmotifs 10 -minw 6 -maxw 21 -objfun de -neg ../5_PauseLoci/${mut}.IDRrep.pauseWindowsSHUFFLE.fasta  -revcomp -markov_order 0

    # Copy results into separate directory
    cp ${mut}/meme.html MEMEout/${mut}_meme.html

    # Get list of significant motifs in separate file
    cat ${mut}/meme.txt | grep MOTIF | grep E-value | sed -E -e 's/[[:blank:]]+/\t/g' | awk -F "\t" -v mut="$mut" '($18 < 0.05) {print mut, $2, $3, mut"_"$3}' | sed 's/ /\t/g' > ${mut}.tmpMotifs.txt

    # Run TOMTOM to find matching motifs
    ./Scripts/matchMotifs.sh ${mut}.tmpMotifs.txt
    
done < $1

