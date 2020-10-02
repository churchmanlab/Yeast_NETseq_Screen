#!/usr/bin/env python

# This scirpt will remove reads that could be generated from splice intermediates that map the either the
# 3' or 5' SS of an intron

# Written by: K. Lachance (inspired by script by H.L.)

# Date: November 8, 2018

# Use: python removeSIfromBAM.py iBAM iSI oBAM

#import pysam

import sys, pysam, os, numpy, re

iBAM = pysam.Samfile(sys.argv[1], 'rb')
iSI = set(open(sys.argv[2], 'r').readlines())
oBAM = pysam.Samfile(sys.argv[3], 'wb', template=iBAM)

MB = set()

# read through starting bam file
for read in iBAM:
    chrom = iBAM.getrname(read.tid)
    
    # selecting the 3' and 5' positions for pos strand 
    if read.is_reverse:
        start = read.aend
        std='pos'
        coord3 = str(chrom)+"_"+str(std)+"_"+str(start)+"\n"

        end = read.pos
        coord5 = str(chrom)+"_"+str(std)+"_"+str(end+1)+"\n"

        
    # selecting the 3' and 5' positions for neg strand 
    if not read.is_reverse:
        start = read.pos
        std='neg'
        coord3 = str(chrom)+"_"+str(std)+"_"+str(start+1)+"\n"

        end = read.aend
        coord5 = str(chrom)+"_"+str(std)+"_"+str(end)+"\n"
        
    # Write output not containing splicing intermediates
    if (coord3 not in iSI) and (coord5 not in iSI):
        oBAM.write(read)

            
iBAM.close()
oBAM.close()
