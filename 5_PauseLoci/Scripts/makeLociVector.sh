#!/bin/bash

# Calcualte number overlapping pauses between all replicates
# Use: ./makeLociVector.sh deletionStrains.txt

# Loop through list of muts twice

while read mut1
do

    while read mut2
    do

	n=`grep -f ${mut1}.IDRrepPausesAll.bed ${mut2}.IDRrepPausesAll.bed | wc -l`
	N=`cat ${mut1}.IDRrepPausesAll.bed | wc -l`
	p=`awk -v n=$n -v N=$N 'BEGIN {print (n/N)*100}'`

	echo -e "${mut1}\t${mut2}\t${n}\t${N}\t${p}"  >> pauseLociVector.txt

    done < $1

done < $1
