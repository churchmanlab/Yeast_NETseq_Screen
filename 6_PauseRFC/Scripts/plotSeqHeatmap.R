#!/n/app/R/3.6.1/bin/Rscript

# Plot heatmap of trinucleotide enrichment

# Written by: K. Lachance
# Date: May 11, 2020

# Use: ./plotSeqHeatmap

# Load libraries
cat("Loading libraries...\n")
library(ggplot2)
library(svglite)
library(reshape2)

data = read.table("PauseSeqFreq.txt", stringsAsFactors = FALSE)
colnames(data) = c("Tri", "N1", "N2", "N3", "Value", "Category")

out = data.frame(matrix(NA, nrow = 16, ncol = 4), row.names = c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"), stringsAsFactors = FALSE)
colnames(out) = c("A", "C", "G", "T")

# Loop through trinucleotides
Nucs = c("A", "C", "G", "T")

for (n1 in Nucs) {
  for (n2 in Nucs) {
    for (n3 in Nucs) {
      
      R = paste(n1, n3, sep="")
      C = n2
      TRI = paste(n1, n2, n3, sep="")
      
      out[R, C] = log(data[data$Tri == TRI & data$Category == "Real", 'Value'] / data[data$Tri == TRI & data$Category == "Shuffle", 'Value'], 2)
      
    } 
  }
}

mDF <- melt(as.matrix(out)) # reshape into dataframe
names(mDF)[c(1:2)] <- c("Y", "X")

p = ggplot(data = mDF, aes(x = X, y = Y, fill = value)) + 
  geom_raster() +
  scale_fill_gradient2(low="grey50", high= "navy", mid = "white", midpoint = 0) + 
  scale_y_discrete(limits = rev(levels(mDF$Y))) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


ggsave(filename = "TriNucSeq_Heatmap.svg", plot = p, height = 4, width = 3)