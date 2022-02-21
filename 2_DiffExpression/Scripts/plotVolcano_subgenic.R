#!/n/app/R/4.0.1/bin/Rscript

# Plot volcano plots for specific deletion strains

# Written by: K. Lachance
# Date: March 30, 2020

# Use: plotVolcano.R mut
# muts="bre1 bur2 bye1 cac1 cac2 cac3 cbc1 ccr4 cdc39 cdc73 chd1 ctr9 dhh1 dst1 eaf1 elf1 gcn5 hda1 hpc2 htz1 ino80 isw1 isw2 leo1 nap1 npl3 paf1 rad6 rco1 rpb4 rph1 rsc30 rtf1 rtt103 set2 set3 spt4 swr1 ubp8 vps15 vps34"
# for mut in $muts
# do
# sbatch -p short -t 0-01:00 -o logs/volcano_sub_${mut}_noAnti.log -e logs/volcano_sub_${mut}_noAnti.err --wrap="../Scripts/plotVolcano_subgenic.R $mut"
# done
# 
# Load libraries
library(ggplot2)
library(reshape2)
library(gridExtra)

# Read in parameters
args <- commandArgs(trailingOnly = TRUE)
mut <- args[1]

# Read in data
data <- read.table(paste(mut, "_subgenic.ALLgenesNoAnti.txt", sep=""), header = TRUE, row.names = 1, stringsAsFactors = FALSE) # ".ALLgenes.txt"

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

ggsave(file=paste(mut, "_subgenic_volcanoNoAnti.pdf", sep=""), plot=p1, width=5, height=5)
