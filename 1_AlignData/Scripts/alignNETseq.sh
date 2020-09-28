#!/bin/bash

# Script to align NET-seq reads to the yeast genome

# Optimized to run on HMS O2 computing cluster (SLURM job scheduler)

# Written by: K. Lachance (modified from similiar scripts by K.H. and B.S.)
# Date: October 29, 2018

# Use: ./alignNETseq.sh alignmentParameters.txt
# NOTE: Run this script a total of 6 times, uncommenting a single block for each run to avoid dependency problems as jobs run
# Leave everything above BLOCK 1 uncommented for each run

# Requires input file with all necessary alignment information (see sample for file format)
parameterFile=$1
project=`grep "Project Name" $parameterFile | perl -pe 's/^.*?(:)//' | sed -e 's/^[ \t]*//' -e 's/[ \t]*$//'`
sampleList=`grep "Sample Names" $parameterFile | perl -pe 's/^.*?(:)//' | sed -e 's/^[ \t]*//' -e 's/[ \t]*$//'`
samples=`echo -e $sampleList`
email=`grep "Notification Email" $parameterFile | perl -pe 's/^.*?(:)//' | sed -e 's/^[ \t]*//' -e 's/[ \t]*$//'`
SIfile=`grep "Splicing Intermediate file" $parameterFile | perl -pe 's/^.*?(:)//' | sed -e 's/^[ \t]*//' -e 's/[ \t]*$//'`
index=`grep "Index Directory" $parameterFile | perl -pe 's/^.*?(:)//' | sed -e 's/^[ \t]*//' -e 's/[ \t]*$//'`

for s in $samples
do

    echo "Working on sample ${s}"
    
    outDir=./postCleaning/${s}
    mkdir -p $outDir/LogErr
    mkdir -p $outDir/Removed

    inDir=../FastqFiles
    in1=${s}.fastq # Files should already be unzipped

    preout1=${outDir}/${s}_noAdaptR_1.fastq
    good1=${s}_cleaned_1
    bad1=${outDir}/Removed/${s}_removed_1

    Adapter=`less ./cutadaptFiles.txt`
    prinseq=./Scripts/prinseq-lite.pl

    ### BLOCK 1: preparing reads with cutadapt, prinseq, and removing UMI
    # Note that files ending in '_noAdaptR_1' are from cutadapt output, '_cleaned_1' are from prinseq, 
    # and '_cleaned_1_noBarcode' are from extract molecular barcode

    # Cut adapter
    sbatch -t 0-4:00 -J ${s}_cutadapt1 -o ${outDir}/LogErr/${s}_cutadapt1.log -e ${outDir}/LogErr/${s}_cutadapt1.err -p short --mail-user=$email --mail-type=ALL cutadapt -f fastq -a $Adapter -O 3 -m 1 --length-tag 'length=' -o ${preout1} ${inDir}/${in1}

    # Clean-up from previous step
    sbatch --mem-per-cpu=5G -J ${s}_cleaned1 -o ${outDir}/LogErr/${s}_cleaned1.log -e ${outDir}/LogErr/${s}_cleaned1.%a.err -t 0-4:00 -p short --mail-user=$email --mail-type=ALL ${prinseq} -fastq ${preout1} -out_good ${outDir}/${good1} -out_bad ${bad1} -no_qual_header -min_len 7 -min_qual_mean 20 -trim_right 1 -trim_ns_right 1 -trim_qual_right 20 -trim_qual_type mean -trim_qual_window 5 -trim_qual_step 1
 
    # Extract barcode
    sbatch -J ${s}_extractBarcode -o ${outDir}/LogErr/${s}_extractBarcode.log -e ${outDir}/LogErr/${s}_extractBarcode.err -t 0-1:00 -p short --mail-user=$email --mail-type=ALL ./Scripts/extractMolecularBarcode.py ${outDir}/${good1}.fastq ${outDir}/${good1}_noBarcode.fastq ${outDir}/${s}_barcodeDistribution.txt ${outDir}/${s}_ligationDistribution.txt

    ### BLOCK 2: running Tophat with reads WITHOUT UMI
    # Note that align_summary.txt will contain alignment statistics

    outDir=./TopHat2/wRTbias/${s}/GENCODE
    mkdir -p $outDir/LogErr
    key=TH_${s}_GENCODE
    inDir=./postCleaning
    reads=${inDir}/${s}/${s}_cleaned_1_noBarcode.fastq

    seg=20
    mis=3

   sbatch --mem-per-cpu=5G -J ${key} -o ${outDir}/LogErr/${key}.log -e ${outDir}/LogErr/${key}.err -t 0-12:00 -p short --mail-user=$email --mail-type=ALL tophat --read-mismatches ${mis} --read-gap-length 2 --read-edit-dist ${mis} -o $outDir --min-anchor-length 8 --splice-mismatches 1 --min-intron-length 50 --max-intron-length 1200 --max-insertion-length 3 --max-deletion-length 3 --num-threads 4 --max-multihits 100 --library-type fr-firststrand --segment-mismatches ${mis} --no-coverage-search --segment-length ${seg} --min-coverage-intron 50 --max-coverage-intron 100000 --min-segment-intron 50 --max-segment-intron 500000 --b2-sensitive $index ${reads}

    ### BLOCK 3: running Tophat with reads WITH UMI

    outDir=./TopHat2/wRTbias/${s}/GENCODE_withBC
    mkdir -p $outDir/LogErr
    reads=${inDir}/${s}/${s}_cleaned_1.fastq

    sbatch --mem-per-cpu=5G -J ${key}_withBC -o ${outDir}/LogErr/${key}.log -e ${outDir}/LogErr/${key}.err -t 0-12:00 -p short --mail-user=$email --mail-type=ALL tophat --read-mismatches ${mis} --read-gap-length 2 --read-edit-dist ${mis} -o $outDir --min-anchor-length 8 --splice-mismatches 1 --min-intron-length 50 --max-intron-length 1200 --max-insertion-length 3 --max-deletion-length 3 --num-threads 4 --max-multihits 100 --library-type fr-firststrand --segment-mismatches ${mis} --no-coverage-search --segment-length ${seg} --min-coverage-intron 50 --max-coverage-intron 100000 --min-segment-intron 50 --max-segment-intron 500000 --b2-sensitive $index ${reads}

    ### BLOCK 4: removing RT mispriming / RT artefact reads
    # Note that useful reads are in accepted_hits_noRTbias_sorted.bam
    
    script=./Scripts/removeReadsWithMismatchesIn6FirstNct.py
    comDir=./TopHat2/wRTbias/
    unmapped='unmapped/'

    samtools view -H ${comDir}/${s}/GENCODE/accepted_hits.bam > ${comDir}/headers.sam

    sbatch -p short -t 0-12:00 -J TH${s}_noBias -o ${outDir}/LogErr/${s}_noBias.log -e ${outDir}/LogErr/${s}_noBias.err --mail-user=$email --mail-type=ALL --wrap="samtools view ${comDir}/${s}/GENCODE/accepted_hits.bam | sed 's/_MolecularBarcode/\t_MolecularBarcode/' | sort -k1,1  > ${comDir}/${s}/GENCODE/accepted_hits_readNamesSorted.sam ; samtools view -h ${comDir}/${s}/GENCODE_withBC/accepted_hits.bam | python ${script} | cut -f 1 | sort -k1,1 | uniq > ${comDir}/${s}/GENCODE/rtBias_readNamesSorted.txt; join -v 1 ${comDir}/${s}/GENCODE/accepted_hits_readNamesSorted.sam ${comDir}/${s}/GENCODE/rtBias_readNamesSorted.txt | sed -e 's/ _MolecularBarcode/_MolecularBarcode/' -e 's/ /\t/g' | cat ${comDir}/headers.sam - | samtools view -Sb -o ${comDir}/${s}/GENCODE/accepted_hits_noRTbias.bam - ; samtools sort ${comDir}/${s}/GENCODE/accepted_hits_noRTbias.bam -o ${comDir}/${s}/GENCODE/accepted_hits_noRTbias_sorted.bam ; samtools index ${comDir}/${s}/GENCODE/accepted_hits_noRTbias_sorted.bam ; rm ${comDir}/${s}/GENCODE/accepted_hits_readNamesSorted.sam ; rm ${comDir}/${s}/GENCODE/accepted_hits_noRTbias.bam"

    ### BLOCK 5: removing splicing intermediates
    
    mkdir -p ${comDir}/${s}/GENCODE/noSI/LogErr
    script=./Scripts/removeSIfromBAM.py

    sbatch -p priority -t 0-12:00 --mail-user=$email --mail-type=ALL --mem=10000 --wrap="samtools view ${comDir}/${s}/GENCODE/accepted_hits_noRTbias_sorted.bam | ./Scripts/ncbiToChrSCer.sh - ${comDir}/${s}/GENCODE/tmp1.sam"

    awk 'NR==FNR {a[$1]=$2;next} {for ( i in a) gsub(i,a[i])}1' ./Genomes/Scerevisiae_R64/chrSCer.txt ${comDir}/headers.sam > ${comDir}/headers_chr.sam

    cat ${comDir}/headers_chr.sam ${comDir}/${s}/GENCODE/tmp1.sam | sed 's/[ \t]\+$//' > ${comDir}/${s}/GENCODE/tmp2.sam

    sbatch -p priority -t 0-12:00 --mail-user=$email --mail-type=ALL --mem=500M --wrap="samtools view -t -h -S -b ${comDir}/${s}/GENCODE/tmp2.sam > ${comDir}/${s}/GENCODE/accepted_hits_noRTbias_chr_sorted.bam"

    sbatch -p short -t 0-12:00 -J ${s}_noSI -o ${comDir}/${s}/GENCODE/noSI/LogErr/${s}_noSI.log -e ${comDir}/${s}/GENCODE/noSI/LogErr/${s}_noSI.err --mail-user=$email --mail-type=ALL $script ${comDir}/${s}/GENCODE/accepted_hits_noRTbias_chr_sorted.bam $SIfile ${comDir}/${s}/GENCODE/accepted_hits_noRTbias_noSI_chr_sorted.bam

    ### BLOCK 6: making bedgraphs from accepted reads

    mkdir -p ./Coverage_noBias/customGenome/${s}/LogErr
    inBam=./TopHat2/wRTbias/${s}/GENCODE/accepted_hits_noRTbias_noSI_chr_sorted.bam
    uniqKey=${s}
    script=./Scripts/customCoverage.py

    rm -fr ${comDir}/${s}/GENCODE/tmp1.sam ${comDir}/${s}/GENCODE/tmp2.sam

    sbatch --mem-per-cpu=5G -t 0-5:00 -J ${uniqKey}_UniqCov -o ./Coverage_noBias/customGenome/${s}/LogErr/${uniqKey}_uniqCustomCov.log -e ./Coverage_noBias/customGenome/${s}/LogErr/${uniqKey}_both_uniqCustomCov.err -p short --mail-user=$email --mail-type=ALL --wrap="samtools view -q 50 ${inBam} | python $script ${uniqKey}"

done
