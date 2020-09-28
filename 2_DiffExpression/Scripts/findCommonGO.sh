#!/bin/bash

while read go
do
    go_id=`echo $go | awk -F "GO:" '{print $2}' | sed 's/)//g'`
    mut=`grep $go_id *_GO.txt | awk -F "_" '{print $1}'`
    
    echo -e $go"\t"$mut >> GOoverlaps.txt

done < ALL.GO.txt 
