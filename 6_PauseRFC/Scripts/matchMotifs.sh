#!/bin/bash

module load gcc/6.2.0
module load openmpi/3.1.0
module load meme/5.0.3

# Run TomTom to find overlapping motifs
mkdir -p TOMTOM

# Read through each motif
while read motif
do

    # Get motif information
    mut=`echo $motif | cut -d " " -f 1`
    id=`echo $motif | cut -d " " -f 2`
    nMotif=`echo $motif | cut -d " " -f 3`
    comb=`echo $motif | cut -d " " -f 4`

    tomtom -evalue -thresh 0.1 -m $id -oc ${mut}/${comb} ${mut}/meme.txt ../0_Annotations/motif_databases/YEAST/YEASTRACT_20130918.meme 
    cp ${mut}/${comb}/tomtom.tsv TOMTOM/${comb}.out

done < $1
