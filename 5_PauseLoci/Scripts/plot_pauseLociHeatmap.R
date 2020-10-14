#!/n/app/R/3.6.1/bin/Rscript

# Plot clustered heatmap with most percent pauses overlapping

# Written by: K. Lachance
# Date: June 4, 2019

# Install and load libraries
cat("Loading libraries...\n")
library("ggplot2")
library("reshape2")
library("hybridHclust")

# Read files
df = read.table("pauseLociVector.txt", row.names = 1, stringsAsFactors = FALSE)

# Only use positions with pauses in at least 2 deletion strains
countOne = apply(df, 2, function(x) sum(x == 1) > 2)
df = df[,countOne]

# Remove poistions where there are no -1s (signal, no pause)
suffSig = apply(df, 2, function(x) sum(x == -1) > 0)
df = df[, suffSig]

# Convert 0 (no data) to NA
df[df == 0] = NA

cat("Plotting heatmap...\n")

# Plot heatmap

# Cluster
row.order <- eisenCluster(df, method = "correlation", compatible = FALSE, verbose = FALSE)$order
col.order <- eisenCluster(t(df), method = "correlation", compatible = FALSE, verbose = FALSE)$order

df2 <- df[row.order, col.order]
m_df2 <- melt(as.matrix(df2))
names(m_df2)[c(1:2)] <- c("Mut", "Pause")

# Plot
p = ggplot(data = m_df2, aes(x = Mut, y = Pause, fill = value)) +
 geom_tile() +
scale_fill_gradient2(name=NULL, low="white", high="darkcyan", na.value = "grey50", labels = NULL) +
  theme_bw() +
theme(axis.title=element_blank(), axis.text.y=element_blank(), axis.ticks.y = element_blank()) +
theme(axis.text.x = element_text(angle = 90)) + 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
   theme(legend.position="none")

ggsave(file="pauseHeatmap.pdf", plot=p, width=8, height=10)
