#!/n/app/R/3.6.1/bin/Rscript

# Make a heatmap for a gene x mutant AS/S ratio matrix

# Written by: K. Lachance
# Date: October 30, 2018

# Use: ./plotWholeGeneHeatmap.R mut

# Load libraries
library(ggplot2)
library(reshape2)
library(svglite)
library(zoo)

# Get arguments (matrices)
args <- commandArgs(trailingOnly = TRUE)
mut <- args[1]

# Load data
cat("Loading data...\n")
mut_pos = read.table(paste(mut, "_OverlapMatrix.WHOLE4250.pos.txt", sep=""), stringsAsFactors = FALSE)
mut_neg = read.table(paste(mut, "_OverlapMatrix.WHOLE4250.neg.txt", sep=""), stringsAsFactors =	FALSE)
mat_mut = rbind(mut_pos, mut_neg)

wt_pos = read.table("wt_OverlapMatrix.WHOLE4250.pos.txt", stringsAsFactors = FALSE)
wt_neg = read.table("wt_OverlapMatrix.WHOLE4250.neg.txt", , stringsAsFactors = FALSE)
mat_wt = rbind(wt_pos, wt_neg)
       
# Save gene information for later
lengths = mat_wt$V3 - mat_wt$V2

# Add pseudocount
cat("Normalizing...\n")
mat_wt <- mat_wt[,5:ncol(mat_wt)] + 0.001
mat_mut <- mat_mut[,5:ncol(mat_mut)] + 0.001

# Normalize
norm_mut = log2(mat_mut / mat_wt)

# Sort by gene size
cat("Sorting by gene length...\n")
order_df = rev(order(lengths))
mat = norm_mut[order_df,]

# Roll apply mean 
cat("Smoothing...\n")
smooth_mat = rollapply(data = t(mat), width = 10, FUN = mean, by.column = TRUE, fill = NA, align = "left")
mat = t(smooth_mat)

# Melt
cat("Melting data frame...\n")
colnames(mat) = 1:ncol(mat)
row.names(mat) = 1:nrow(mat)
melted_mat = melt(as.matrix(mat))
colnames(melted_mat) = c("Gene", "Position", "Value")

# Plot
cat("Plotting...\n")
p <- ggplot(melted_mat, aes(x = Position, y = Gene, fill = Value)) + 
  geom_tile() + 
  scale_fill_gradient2(name=NULL, limits = c(-11, 15), low="#A7533F", mid="white", high="#507BA7", midpoint=0, na.value = "white") + 
  geom_vline(xintercept = 250, linetype="dotted") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank()) + 
  theme(panel.background = element_blank()) + 
  theme(legend.position = "none")

ggsave(paste(mut, "_wholeGeneASHeatmap.pdf", sep=""), p, width = 2, height = 3)