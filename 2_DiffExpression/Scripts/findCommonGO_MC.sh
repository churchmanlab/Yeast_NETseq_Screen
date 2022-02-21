#!/bin/bash

# Use ../Scripts/findCommonGO_MC.sh
# sbatch -p short -t 0-01:00 --wrap="./Scripts/findCommonGO_MC.sh"


rm GOoverlapsUP.txt
rm GOoverlapsDOWN.txt
rm GOoverlapsALL.txt

while read goup
do
    go_id=`echo $goup`
    mut=`grep "^${go_id}$" *_UP_GOterms.txt | awk -F "_" '{print $1}'`
    
    echo -e $goup"\t"$mut >> GOoverlapsUP.txt

done < ALL.UP_GO.txt 

while read godown
do
    go_id=`echo $godown`
    mut=`grep "^${go_id}$" *_DOWN_GOterms.txt | awk -F "_" '{print $1}'`
    
    echo -e $godown"\t"$mut >> GOoverlapsDOWN.txt

done < ALL.DOWN_GO.txt 

while read goall
do
    go_id=`echo $goall`
    mut=`grep "^${go_id}$" *_ALL_GOterms.txt | awk -F "_" '{print $1}'`
    
    echo -e $goall"\t"$mut >> GOoverlapsALL.txt

done < ALL.ALL_GO.txt 

sed -i s/\'//g GOoverlaps*
