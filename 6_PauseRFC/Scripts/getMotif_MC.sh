#!/bin/bash

#SBATCH -c 4
#SBATCH -t 0-01:00
#SBATCH -p short
#SBATCH --mem=80G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=youremail@gmail.com

###
###

# Use: 
# sbatch -o logs/6_getMotif.log -e logs/6_getMotif.err ./Scripts/getMotif_MC.sh deletionStrainsCovFiltered.txt



# Must load all in this order to work
module purge
module load gcc/6.2.0
module load ghostscript/9.24
echo loading openmpi... >> logs/6_getMotif.err
module load openmpi/3.1.0
echo done loading openmpi >> logs/6_getMotif.err
echo loading meme... >> logs/6_getMotif.err
module load meme/5.0.3
echo done loading meme >> logs/6_getMotif.err

mkdir MEMEout

while read mut
do

    echo $mut
    mkdir -p ${mut}

	# Shuffle pause windows for control file
# 	echo shuffling pause windows for ${mut}... > logs/6_getMotif.err
# 	fasta-shuffle-letters ../5_PauseLoci/${mut}.IDRrep.pauseWindows.fasta ../5_PauseLoci/${mut}.IDRrep.pauseWindowsSHUFFLE.fasta
# 	echo done shuffling pause windows for ${mut} > logs/6_getMotif.err
	
    # Get motifs
    # NOTE: must not be on interactive mode for this to work
    echo running meme for ${mut}... >> logs/6_getMotif.err
    meme ../5_PauseLoci/${mut}.IDRrep.pauseWindows.fasta -dna -oc $mut -time 18000 -nostatus -mod zoops -nmotifs 10 -minw 6 -maxw 21 -objfun de -neg ../5_PauseLoci/${mut}.IDRrep.pauseWindowsSHUFFLE.fasta -revcomp -markov_order 0 # --with-gs=/n/app/ghostscript/9.24/bin/gs
	echo done running meme for ${mut}... >> logs/6_getMotif.err
	
	
    # Copy results into separate directory
    cp ${mut}/meme.html MEMEout/${mut}_meme.html

    # Get list of significant motifs in separate file
    cat ${mut}/meme.txt | grep MOTIF | grep E-value | sed -E -e 's/[[:blank:]]+/\t/g' | awk -F "\t" -v mut="$mut" '($18 < 0.05) {print mut, $2, $3, mut"_"$3}' | sed 's/ /\t/g' > ${mut}.tmpMotifs.txt

    # Run TOMTOM to find matching motifs
    ./Scripts/matchMotifs.sh ${mut}.tmpMotifs.txt
    
done < $1

