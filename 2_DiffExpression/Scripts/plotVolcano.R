#!/n/app/R/3.6.1/bin/Rscript

# Plot volcano plots for specific deletion strains

# Written by: K. Lachance
# Date: March 30, 2020

# Use: plotVolcano.R mut

# Load libraries
library(ggplot2)
library(reshape2)
library(gridExtra)
library(svglite)

# Read in parameters
args <- commandArgs(trailingOnly = TRUE)
mut <- args[1]

# Read in data
data <- read.table(paste(mut, ".ALLgenes.txt", sep=""), header = TRUE, row.names = 1, stringsAsFactors = FALSE)

# Remove rows with NA
data = data[!is.na(data$log2FoldChange) & !is.na(data$padj),]

# Add color code
data$Col = "G"
data[data$log2FoldChange > 1 & data$padj < 0.05, 'Col'] = "B"
data[data$log2FoldChange < -1 & data$padj < 0.05, 'Col'] = "R"

nUp = nrow(data[data$Col == "B",])
nDown = nrow(data[data$Col == "R",])

cat(paste(mut, " has ", nUp, " genes up-regulated and ", nDown, " genes down-regulated, for a total of ", nUp + nDown, " total differentially expressed genes!\n", sep = ""))

# Plot
p1 <- ggplot(data, aes(x = log2FoldChange, y = -log10(padj), color = Col, fill = Col)) + 
  geom_point(size = 1) + 
  scale_color_manual(values = c("G" = "grey50", "B" = "dodgerblue3", "R" = "coral3")) + 
  scale_fill_manual(values = c("G" = "grey50", "B" = "dodgerblue3", "R" = "coral3")) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(legend.position = "none")

ggsave(file=paste(mut, "_volcano.svg", sep=""), plot=p1, width=5, height=5)
