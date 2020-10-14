# Yeast NETseq Screen
Dynamics of transcription elongation are finely-tuned by dozens of regulatory factors

![Screen Overview](https://github.com/churchmanlab/Yeast_NETseq_Screen/blob/master/ScreenOverview.png)

This repository includes the scripts and annotation files needed to analyze NET-seq data generated from *S. cerevisiae*. The directories are listed in order to take raw NET-seq Fastq files through alignment and all analyses presented in the paper in order. There are README files for each analysis step explaining how scripts are run and what is generated with each. 

# Analysis steps
1. Align NET-seq data
2. Identify differentially expressed genes and enriched gene ontology terms
3. Quantify antisense transcription
4. Calculate pausing indices around RNA processing sites
5. Identify Pol II pause loci
6. Build a random forest classifier to predict Pol II pause loci
7. Correlated transcriptional phenotypes with one another

# Data availability and manuscript
Fastq files for NET-seq data will be deposited to GEO; when available, the accession number will be listed here. The link to our full manuscript will be provided here upon publication. 
