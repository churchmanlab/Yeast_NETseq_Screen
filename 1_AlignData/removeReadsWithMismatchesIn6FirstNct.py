#!/usr/bin/env python

"""

Date : May 9th 2014

Script used to filter reads that contain mismatches in the 6 first nucleotides (that corresponds to the barcode if it
hasn't been removed)

use : python removeReadsWithMismatchesIn6FirstNct.py stdin stdout

"""

import sys, pysam, re


iBAM = pysam.Samfile("-", 'r')
oBAM = pysam.Samfile("-", 'w', template=iBAM)

for line in iBAM:
    md=re.findall(r'\d+', [tag[1] for tag in line.tags if tag[0]=='MD'][0])
    if len(md) == 1 :
        oBAM.write(line)
    else:
        if (not line.is_reverse) and (int(md[0]) > 5):
            oBAM.write(line)
        elif (line.is_reverse) and (int(md[-1]) > 5):
            oBAM.write(line)

iBAM.close()
oBAM.close()

