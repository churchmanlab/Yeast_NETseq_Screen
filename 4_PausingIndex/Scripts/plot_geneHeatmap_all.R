#!/n/app/R/3.6.1/bin/Rscript

# Written by: K. Lachance
# Date: July 20, 2020

# Make a heatmap for a gene x mutant PI

# Use: ./plot_geneHeatmap_all.R REGION

# Load libraries
library(ggplot2)
library(reshape2)
library(svglite)

# Get arguments (matrices)
args <- commandArgs(trailingOnly = TRUE)
REGION <- args[1]

# Get data by combining tables
data_all = read.table("ALL.PI.bed", stringsAsFactors = FALSE)
colnames(data_all) = c("Chr", "Start1", "End1", "Start2", "End2", "Gene", "Skip", "Strand", "ROIcount", "notROIcount", "ROInorm", "notROInorm", "PI", "Region", "Mut")

data = data_all[data_all$Region == REGION,]

# Get data in form for clustering
muts = unique(data$Mut)
genes = unique(data$Gene)
out = data.frame(matrix(0, nrow = length(muts), ncol = length(genes)), stringsAsFactors = FALSE)
colnames(out) = genes
row.names(out) = muts
for (m in muts) {
    for (g in genes) {
	PI = data[data$Gene == g & data$Mut == m, 'PI']
	if (length(PI) > 0) { out[m, g] = PI } 
	else { out[m, g] = NA }
    }
}

# Remove columns with NA / Inf
out = out[,is.finite(colSums(out))]

# Cluster
row.order <- hclust(dist(out, method = "euclidean"), method="ward.D")$order
col.order <- hclust(dist(t(out), method = "euclidean"), method="ward.D")$order
out_new <- out[row.order, col.order]
m_out <- melt(as.matrix(out_new))
names(m_out)[c(1:2)] <- c("Mut", "Gene")

# Set color
if (REGION == "TSS") {col = "#006838"}
if (REGION == "pA") {col = "#BE1E2D"}
if (REGION == "5pSS") {col = "#2E3192"}
if (REGION == "3pSS") {col = "#00AEEF"}
if (REGION == "Anti") {col = "#507BA7"}

# Plot
p1 <- ggplot(m_out, aes(x=Gene, y=Mut, fill=log2(value))) +
   geom_tile() +
   scale_fill_gradient2(low = "grey50", mid = "white", high = col) + 
   theme_bw() +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
   theme(axis.text.x = element_text(angle = 90, vjust = 1, size = 4, hjust = 1)) +
   theme(axis.text.y = element_text(size=4)) +
   theme(axis.title=element_blank(), axis.text.x = element_blank())

filename = paste("ALL_", REGION, "_GenePIheatmap.svg", sep="")
cat(paste("Saving plot to: ", filename, "\n", sep=""))

ggsave(file=filename, plot=p1, width=7, height=6)
