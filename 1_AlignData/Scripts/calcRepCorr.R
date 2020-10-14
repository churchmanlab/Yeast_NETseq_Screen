#!/n/app/R/3.6.1/bin/Rscript

# Calculate gene-by-gene RPKM correlation between two replicates

# Written by: K. Lachance
# Date: November 26, 2018

# Use: ./calcRepCorr.R mut-1 mut-2

# Install and load libraries
library(ggplot2)
library(svglite)

# Get input parameters
args <- commandArgs(trailingOnly = TRUE)
mut1_name <- args[1]
mut2_name <- args[2]

# Get full file name from mutant names
mut1_file = paste("./RPKMnormFiles/", mut1_name, ".RPKMnorm.bed", sep="")
mut2_file = paste("./RPKMnormFiles", mut2_name, ".RPKMnorm.bed", sep="")

# Read tables
mut1 <- read.table(mut1_file)
mut2 <- read.table(mut2_file)

# Combine tables
dat <- cbind(mut1, mut2)[,c(1:7, 14)]
colnames(dat) <- c("Chr", "Start", "End", "Gene", "Label", "Strand", "RPKM1", "RPKM2")

# Calculate RPKM correlation
repCor=as.character(round(cor(dat$RPKM1, dat$RPKM2, method='pearson'), digits=2))
cat(paste("Gene-by-Gene RPKM Correlation between ", mut1_name, " and ", mut2_name, ": ", repCor, sep=""))

# Plot
filename = paste("repCorr.", mut1_name, ".", mut2_name, ".svg", sep="")

p <- ggplot(dat, aes(x=log10(RPKM1), y=log10(RPKM2))) + 
  geom_point(alpha = 0.4) + 
  geom_abline(intercept=0, slope=1, col='red', linetype=2) + 
  theme_bw() + 
  xlab(paste("log10(", mut1_name, " RPKM)", sep="")) + 
  ylab(paste("log10(", mut2_name, " RPKM)", sep="")) + 
  xlim(-2, 4) + ylim(-2, 4) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(file=filename, plot=p, width=6, height=6)
