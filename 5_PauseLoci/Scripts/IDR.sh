#!/bin/bash

# Find pauses above IDR threshold for each mut

while read mut
do

    echo $mut

    # Get list of pauses with scores
    ./Scripts/findOverlap_pauseScore.R $mut

    # Calculate and plot IDR
    ./Scripts/plot_pauseScore.R $mut


#done < deletionStrains.txt
done < tmp
