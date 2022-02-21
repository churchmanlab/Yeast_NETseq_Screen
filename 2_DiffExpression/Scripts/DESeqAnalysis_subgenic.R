#!/n/app/R/4.0.1/bin/Rscript

# Identify differentially expressed genes

# Written by: K. Lachance
# Date: November 27, 2018
# Updated 12/2021 by M. Couvillion to compare treated vs untreated instead of vice versa

# Use: sbatch -p short -t 0-01:00 --wrap="./Scripts/DESeqAnalysis_subgenic.R strain"

# muts="bre1 bur2 bye1 cac1 cac2 cac3 cbc1 ccr4 cdc39 cdc73 chd1 ctr9 dhh1 dst1 eaf1 elf1 gcn5 hda1 hpc2 htz1 ino80 isw1 isw2 leo1 nap1 npl3 paf1 rad6 rco1 rpb4 rph1 rsc30 rtf1 rtt103 set2 set3 spt4 swr1 ubp8 vps15 vps34"
# for mut in $muts
# do
# sbatch -p short -t 0-01:00 --wrap="../Scripts/DESeqAnalysis_subgenic.R $mut"
# done

# Load libraries
library(DESeq2)
library(apeglm)

# Read in parameters (data)
args <- commandArgs(trailingOnly = TRUE)
strain <- args[1]

# Get data
dat_filename=paste(strain, "_subgenic.mat.txt", sep="")
cond_filename=paste(strain, "_subgenic.cond.txt", sep="")

# Import data as matrices
cts <- as.matrix(read.csv(dat_filename,sep="\t",header=TRUE, row.names="Gene"))

coldata <- read.csv(cond_filename, sep="\t", row.names=1)

# Make sure that columns are named consistantly
colnames(cts) <- row.names(coldata)

# Remove rows with NA from data
cts = cts[complete.cases(cts),]

# Relevel so the comparison is treated to untreated *MTC added 12/2021*
coldata$Condition <- factor(coldata$Condition, levels=c("untreated","treated"))

# Perform differential expression analysis
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Condition)
dds <- DESeq(dds)

# Compute results with alpha = 0.05
res <- results(dds, alpha = 0.05)
resLFC <- lfcShrink(dds, coef="Condition_treated_vs_untreated", type="apeglm")

# Get size factors
# dds <- estimateSizeFactors(dds)
sf <- sizeFactors(dds)
# Get normalized counts
normcts <- counts(dds, normalized=TRUE)

write.table(as.data.frame(res), file=paste(strain, "_subgenic.ALLgenes.txt", sep=""),quote=F, sep="\t")
write.table(as.data.frame(sf), file=paste(strain, "_subgenic.sizeFactors.txt", sep=""),quote=F, col.names=F, sep="\t")
write.table(as.data.frame(normcts), file=paste(strain, "_subgenic.normCounts.txt", sep=""),quote=F, sep="\t")

# Order based on adjusted p-value
resOrdered <- res[order(res$pvalue),]

# Calculate number of DE genes
nDE = sum((res$padj < 0.05 & abs(res$log2FoldChange) > 1), na.rm=TRUE)
cat(paste("**********Number of DE genes for ", strain, ": ", nDE, "**********\n", sep=""))

# Write file with only significantly differentlially expressed genes
resSig <- subset(resOrdered, padj < 0.05)
difSig <- resSig[abs(resSig$log2FoldChange) > 1, ]
resOrdered <- na.omit(resOrdered)
same <- resOrdered[abs(resOrdered$log2FoldChange) < 0.2, ]


write.table(as.data.frame(difSig), file=paste(strain, "_subgenic.DEgenes.txt", sep=""),quote=F, sep="\t")
# write.table(as.data.frame(sameSig), file=paste(strain, ".SAMEgenes.txt", sep=""),quote=F, sep="\t")