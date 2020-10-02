#!/n/app/R/3.6.1/bin/Rscript

# Written by: K. Lachance
# June 8, 2020

# Make a heatmap for a gene x mutant PI

# Example: ./plot_geneHeatmap_geneEnds.R mut

# Load libraries
library(ggplot2)
library(reshape2)
library(svglite)

# Get arguments (matrices)
args <- commandArgs(trailingOnly = TRUE)
mut <- args[1]

# Get data by combining tables
data1 = read.table(paste(mut, "_TSS.PI.bed", sep=""), stringsAsFactors = FALSE)
colnames(data1) = c("Chr", "Start1", "End1", "Start2", "End2", "Gene", "Skip", "Strand", "ROIcount", "notROIcount", "ROInorm", "notROInorm", "PI", "Region", "Mut")
data1 = data1[,c('Gene', 'PI', 'Region')]

data2 =	read.table(paste(mut, "_pA.PI.bed", sep=""), stringsAsFactors = FALSE)
colnames(data2) = c("Chr", "Start1", "End1", "Start2", "End2", "Gene", "Skip", "Strand", "ROIcount", "notROIcount", "ROInorm", "notROInorm", "PI", "Region", "Mut")
data2 =	data2[,c('Gene', 'PI', 'Region')]

data3 = read.table(paste(mut, "_Antisense.PI.bed", sep=""), stringsAsFactors = FALSE)
colnames(data3) = c("Chr", "Start1", "End1", "Start2", "End2", "Gene", "Skip", "Strand", "ROIcount", "notROIcount", "ROInorm", "notROInorm", "PI", "Region", "Mut")
data3 = data3[,c('Gene', 'PI', 'Region')]

data = rbind(data1, data2, data3)

# Get data in form for clustering
genes = unique(data$Gene)
out = data.frame(matrix(0, nrow = 3, ncol = length(genes)), stringsAsFactors = FALSE)
colnames(out) = genes
regions = c("TSS", "pA", "Anti")
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


col1 = "#006838"
col2 = "#BE1E2D"
col3 = "#507BA6"

# Plot
p1 <- ggplot(m_out[m_out$Region == "TSS",], aes(x=Gene, y=Region, fill=log2(value))) +
   geom_tile() +
   scale_fill_gradient2(low = "grey50", mid = "white", high = col1) + 
   theme_bw() +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
   theme(axis.title=element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) +
   theme(panel.background = element_blank())

filename = paste(mut, "_TSS_GenePIheatmap.svg", sep="")
cat(paste("Saving plot to: ", filename, "\n", sep=""))
ggsave(file=filename, plot=p1, width=7, height=1)

p2 <- ggplot(m_out[m_out$Region == "pA",], aes(x=Gene, y=Region, fill=log2(value))) +
   geom_tile() +
   scale_fill_gradient2(low = "grey50", mid = "white", high = col2) +
   theme_bw() +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
   theme(axis.title=element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) +
   theme(panel.background = element_blank())

filename =  paste(mut, "_pA_GenePIheatmap.svg", sep="")
cat(paste("Saving plot to: ", filename, "\n", sep=""))
ggsave(file=filename, plot=p2, width=7, height=1)

p3 <- ggplot(m_out[m_out$Region == "Anti",], aes(x=Gene, y=Region, fill=log2(value))) +
   geom_tile() +
   scale_fill_gradient2(low = "grey50", mid = "white", high = col3) +
   theme_bw() +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
   theme(axis.title=element_blank(), axis.text = element_blank(), axis.ticks = element_blank()) +
   theme(panel.background = element_blank())

filename =  paste(mut, "_Antisense_GenePIheatmap.svg", sep="")
cat(paste("Saving plot to: ", filename, "\n", sep=""))
ggsave(file=filename, plot=p3, width=7, height=1)