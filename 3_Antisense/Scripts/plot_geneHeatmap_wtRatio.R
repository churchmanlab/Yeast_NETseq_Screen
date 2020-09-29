#!/n/app/R/3.6.1/bin/Rscript

# Make a heatmap for a gene x mutant PI

# Written by: K. Lachance
# Date: May 8, 2019

# Use: ./plot_geneHeatmap.R strainOfInterest.txt VALUE

# Load libraries
library(ggplot2)
library(reshape2)
library(svglite)
library("hybridHclust")

# Get data by combining tables
data = read.table("./CodingAntisense/ALL.geneNETseqRatio_withWTratio.bed", stringsAsFactors = FALSE)
colnames(data) = c("Chr", "Start", "End", "Gene", "Skip", "Strand", "Sense", "Anti", "Ratio", "Mut", "WTsense", "WTantisense")

muts = read.table("deletionStrains.txt", stringsAsFactors = FALSE)$V1

# Calc ratio of ratios
data$RatioRatio = data$WTantisense / data$WTsense
c = 13 # Column of ratio of ratios

# Get data in form for clustering
genes = unique(data$Gene)
out = data.frame(matrix(0, nrow = length(muts), ncol = length(genes)), stringsAsFactors = FALSE)
colnames(out) = genes
row.names(out) = muts
for (m in muts) {
    for (g in genes) {
	V = data[data$Gene == g & data$Mut == m, c]
	if (length(V) > 0) { out[m, g] = V } 
	else { out[m, g] = NA }
    }
}

# Replace Inf with NA
out = do.call(data.frame,lapply(out, function(x) replace(x, is.infinite(x),NA)))
row.names(out) = muts

# Replace NA with 0
out[is.na(out)] <- 0

# Only use positions with signal for at least 1 strain
countOne = apply(out, 2, function(x) sum(x) > 0)
out = out[,countOne]

# Cluster normally
row.order <- hclust(dist(out, method = "euclidean"), method="ward.D")$order
col.order <- hclust(dist(t(out), method = "euclidean"), method="ward.D")$order

out_new <- out[row.order, col.order]

m_out <- melt(as.matrix(out_new))
names(m_out)[c(1:2)] <- c("Mut", "Gene")

# Plot
p1 = ggplot(m_out, aes(x=Gene, y=Mut, fill=log2(value))) +
   geom_tile() +
   scale_fill_gradient2(low = "#a7543f", mid = "white", high = "#507ba6", na.value = "white") +
   theme_bw() +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
   theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 4, hjust = 1)) +
   theme(axis.text.y = element_text(size=4))

filename = "ALL_GeneHeatmap_wtRatio.svg"
cat(paste("Saving plot to: ", filename, "\n", sep=""))

ggsave(file=filename, plot=p1, width=7, height=6)
