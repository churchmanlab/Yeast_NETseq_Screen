#!/n/app/R/3.6.1/bin/Rscript

# Identify differentially expressed genes

# Written by: K. Lachance
# Date: November 27, 2018

# Use: ./DESeqAnalysis.R strain

# Load libraries
library(DESeq2)
library(apeglm)

# Read in parameters (data)
args <- commandArgs(trailingOnly = TRUE)
strain <- args[1]

# Get data
dat_filename=paste(strain, ".mat.txt", sep="")
cond_filename=paste(strain, ".cond.txt", sep="")

# Import data as matrices
cts <- as.matrix(read.csv(dat_filename,sep="\t",row.names="Gene"))

coldata <- read.csv(cond_filename, sep="\t", row.names=1)

# Make sure that columns are named consistantly
colnames(cts) <- row.names(coldata)

# Remove rows with NA from data
cts = cts[complete.cases(cts),]

# Perform differential expression analysis
dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Condition)
dds <- DESeq(dds)

# Compute results with alpha = 0.05
res <- results(dds, alpha = 0.05)
resLFC <- lfcShrink(dds, coef="Condition_untreated_vs_treated", type="apeglm")

write.table(as.data.frame(res), file=paste(strain, ".ALLgenes.txt", sep=""),quote=F, sep="\t")

# Order based on adjusted p-value
resOrdered <- res[order(res$pvalue),]

# Calculate number of DE genes
nDE = sum((res$padj < 0.05 & abs(res$log2FoldChange) > 1), na.rm=TRUE)
cat(paste("**********Number of DE genes for ", strain, ": ", nDE, "**********\n", sep=""))

# Write file with only significantly differentlially expressed genes
resSig <- subset(resOrdered, padj < 0.05)
resSig <- resSig[abs(resSig$log2FoldChange) > 1, ]
write.table(as.data.frame(resSig), file=paste(strain, ".DEgenes.txt", sep=""),quote=F, sep="\t")