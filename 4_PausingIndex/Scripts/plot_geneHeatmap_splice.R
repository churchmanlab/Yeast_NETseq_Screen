#!/n/app/R/3.6.1/bin/Rscript

# Written by: K. Lachance
# Date: June 9, 2020

# Make a heatmap for a gene x mutant PI

# Use: ./plot_geneHeatmap_splice.R mut

# Load libraries
library(ggplot2)
library(reshape2)
library(svglite)

# Get arguments (matrices)
args <- commandArgs(trailingOnly = TRUE)
mut <- args[1]

# Get data by combining tables
data1 = read.table(paste(mut, "_5pSS.PI.bed", sep=""), stringsAsFactors = FALSE)
data2 =	read.table(paste(mut, "_3pSS.PI.bed", sep=""), stringsAsFactors = FALSE)
data = rbind(data1, data2)
colnames(data) = c("Chr", "Start1", "End1", "Start2", "End2", "Gene", "Skip", "Strand", "ROIcount", "notROIcount", "ROInorm", "notROInorm", "PI", "Region", "Mut")

# Get data in form for clustering
genes = unique(data$Gene)
out = data.frame(matrix(0, nrow = 2, ncol = length(genes)), stringsAsFactors = FALSE)
colnames(out) = genes
regions = c("5pSS", "3pSS")
row.names(out) = regions
for (r in regions) {
    for (g in genes) {
	PI = data[data$Gene == g & data$Region == r, 'PI']
	if (length(PI) > 0) { out[r, g] = PI } 
	else { out[r, g] = NA }
    }
}

# Remove columns with NA / Inf
out = out[,is.finite(colSums(out))]

# Cluster
col.order <- hclust(dist(t(out), method = "euclidean"), method="ward.D")$order
out_new <- out[, col.order]
m_out <- melt(as.matrix(out_new))
names(m_out)[c(1:2)] <- c("Region", "Gene")

# Set colors
col1 = "#2E3192"
col2 = "#00AEEF"

# Plot
p1 <- ggplot(m_out[m_out$Region == "5pSS",], aes(x=Gene, y=Region, fill=log2(value))) +
   geom_tile() +
   scale_fill_gradient2(low = "grey50", mid = "white", high = col1) + 
   theme_bw() +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
   theme(axis.title=element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) +
   theme(panel.background = element_blank())

filename = paste(mut, "_5pSS_GenePIheatmap.svg", sep="")
cat(paste("Saving plot to: ", filename, "\n", sep=""))
ggsave(file=filename, plot=p1, width=7, height=1)

p2 <- ggplot(m_out[m_out$Region == "3pSS",], aes(x=Gene, y=Region, fill=log2(value))) +
   geom_tile() +
   scale_fill_gradient2(low = "grey50", mid = "white", high = col2) +
   theme_bw() +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
   theme(axis.title=element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) +
   theme(panel.background = element_blank())

filename = paste(mut, "_3pSS_GenePIheatmap.svg", sep="")
cat(paste("Saving plot to: ", filename, "\n", sep=""))
ggsave(file=filename, plot=p2, width=7, height=1)
