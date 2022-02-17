#!/n/app/R/4.0.1/bin/Rscript

# Perform PCA and plot deletion strain pauses 

# Written by: K. Lachance
# Date: July 14, 2019

# Use: ./Scripts/pausePCA.R 

# Install and load libraries
cat("Loading libraries...\n")
library("ggplot2")

# Load data
dat <- read.table("mutPauseVectors.txt", row.names = 1, stringsAsFactors = FALSE)

# PCA
dat_pca = prcomp(t(dat), scale = FALSE)

pcs = data.frame(dat_pca$rotation)

eigs = (dat_pca$sdev)^2
eigs = eigs[-1]
varE = round((eigs / sum(eigs)) * 100, 2)
vE = data.frame(PC = 1:length(varE), varE = varE, stringsAsFactors = FALSE)

q1 <- ggplot(pcs, aes(x = PC1, y = PC2)) + 
   geom_point() + 
   geom_text(aes(label = rownames(pcs)), hjust=0.005, vjust=0.005, size = 2) +
   theme_bw() + 
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
   ggtitle("PC1 v PC2")

ggsave(plot = q1, file = "pausePCA.pdf", device="pdf", width = 5, height = 5)
