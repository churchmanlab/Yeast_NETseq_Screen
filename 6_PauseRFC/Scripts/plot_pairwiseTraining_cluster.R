#!/n/app/R/3.6.1/bin/Rscript

# Plot heatmap showing AUC for transfer learning

# Written by: K. Lachance
# Date: May 8, 2020

# Use: ./plot_paiwiseTraining.R

# Load libraries
cat("Loading libraries...\n")
library(reshape2)
library(ggplot2)
library(svglite)

# Get input parameters
cat("Importing data...\n")
data <- read.table("AUCvalues_pairwise.txt", stringsAsFactors = FALSE)
colnames(data) = c("Train", "Test", "AUC")

# Add pairs
pairs = read.table("AUCvalues.txt", stringsAsFactors = FALSE)
tmp = data.frame(Train = pairs$V2, Test = pairs$V2, AUC = pairs$V3, stringsAsFactors = FALSE)
data = rbind(data, tmp)

# Put into form for clustering
cat("Formatting for clustering...\n")
Muts = sort(unique(data$Test))
df = data.frame(matrix(NA,nrow = length(Muts), ncol = length(Muts)), row.names = Muts, stringsAsFactors = FALSE)
colnames(df) = Muts

for (m1 in Muts) {
    for (m2 in Muts) {
    	df[m1, m2] = data[data$Test == m1 & data$Train == m2, 'AUC']
    }
}

# Cluster
cat("Clustering...\n")
row.order <- hclust(dist(df))$order # clustering
col.order <- hclust(dist(t(df)))$order
dat <- df[row.order, col.order] # re-order matrix accoring to clustering

df_molten_dat <- melt(as.matrix(dat)) # reshape into dataframe
names(df_molten_dat)[c(1:2)] <- c("Test", "Train")

# Plot
cat("Creating heatmap...\n")
p <- ggplot(df_molten_dat, aes(x = Test, y = Train)) +
  geom_tile(aes(fill = value),color = "white") +
  scale_fill_gradient(low = "white", high = "darkorange", na.value = "darkorange4") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

cat("Saving heatmap to pairwiseTrainingHeatmap_cluster.svg...\n")
ggsave(filename = "pairwiseTrainingHeatmap_cluster.svg", plot = p, height = 8, width = 8)
