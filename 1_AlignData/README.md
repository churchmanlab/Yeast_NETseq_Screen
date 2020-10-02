# Processing and alignment of NET-seq data
The adapter sequence (ATCTCGTATGCCGTCTTCTGCTTG) was removed using cutadapt with the following parameters: -O 3 -m 1 --length-tag ‘length=’. Raw fastq files were filtered using PrinSeq (http://prinseq.sourceforge.net/) with the following parameters: -no_qual_header -min_len 7 -min_qual_mean 20 -trim_right 1 -trim_ns_right 1 -trim_qual_right 20 -trim_qual_type mean -trim_qual_window 5 -trim_qual_step 1. Random hexamer linker sequences (the first 6 nucleotides at the 5’ end of the read) were removed using custom Python scripts, but remained associated with the read. Reads were then aligned to the SacCer3 genome obtained from the Saccharomyces Genome Database using the TopHat2 aligner (Kim et al., 2013) with the following parameters: --read-mismatches 3 --read-gap-length 2 --read-edit-dist 3 --min-anchor-length 8 --splice-mismatches 1 --min-intron-length 50 --max-intron-length 1200 --max-insertion-length 3 --max-deletion-length 3 --num-threads --max-multihits 100 --library-type fr-firststrand --segment-mismatches 3 --no-coverage-search --segment-length 20 --min-coverage-intron 50 --max-coverage-intron 100000 --min-segment-intron 50 --max-segment-intron 500000 --b2-sensitive. To avoid any biases toward favoring annotated regions, the alignment was performed without providing a transcriptome. Reverse transcription mispriming events were identified and removed where molecular barcode sequences correspond exactly to the genomic sequence adjacent to the aligned read. With NET-seq, the 5’ end of the sequencing, which corresponds to the 3’ end of the nascent RNA fragment, is recorded with a custom Python script using the HTSeq package (Anders et al., 2015). NET-seq data were normalized by million mapped reads. Replicate correlations were performed comparing RPKM of each gene in each replicates; replicates were considered highly correlated with a Pearson correlation of R ≥ 0.85. Raw NET-seq data of highly correlated replicates were merged, and then re-normalized by million mapped reads. 

# Software requirements for aligning and analyzing NET-seq data (also in requirements.txt)
- gcc (6.2.0)
- python (2.7.12)
- cutadapt (1.14)
- perl (5.24.0)
- fastqc (0.11.5)
- bedops (2.4.30)
- bedtools (2.27.1)
- samtools (1.9)
- bowtie2 (2.2.9)
- tophat (2.1.1)
- htseq (0.9.1)

# To align NET-seq data
1. Inatall and load all required software
2. Place NET-seq Fastq files in FastqFiles directory
3. Edit alignmentParameters.txt. Can change Project Name, Sample Names, and Notification Email
4. Edit `./Scripts/alignNETseq.sh` such that BLOCKs 2-6 are commented out and not run (only one job sumbission is uncommented at a time). Note that this script is writted to be compatible with a SLURM job scheduler. If your cluster cannot interpret these jobs, edit the script to remove the wrapping job submission commands (`sbatch`...)
5. Run the script from the 1_AlignData directory with the command `./Scripts/alignNETseq.sh alignmentParameters.txt`
6. Use the same command more times, commenting in a different job sumbission BLOCK each time, in sequence. When the last job is run, you will have produced the following files for each sample:
   - `SampleName_cov_pos.bedGraph` and `SampleName_cov_neg.bedGraph`, in which the entire read contributes to the coverage
   - `SampleName_stt_pos.bedGraph` and `SampleName_stt_neg.bedGraph`, in which the beginning of the read is mapped (represents Pol II position)
   - `SampleName_end_pos.bedGraph` and `SampleName_end_pos.bedGraph`, in which the end of the read is mapped

# Once the data are aligned
1. Move all completed coverage files to a new directory, `CoverageFiles/` with the command `mv *.bedGraph CoverageFiles` and remove any other temporary output files
2. Combine and normalize for number of replicates with the command `./Scripts/combineReps.sh SampleName`. This will create output in a `CombinedReplicates/` directory. Note that replicates of each sample name are assumed to be in the format `SampleName-1`, `SampleName-2`, etc.
3. Create RPM-normalized files with the command `./Scripts/calcRPM.sh SampleName`. Output files will be in the `RPMnormFiles` directory
4. 

