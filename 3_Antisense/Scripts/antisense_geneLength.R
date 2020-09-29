#!/n/app/R/3.6.1/bin/Rscript

# Written by: K. Lachance
# Date: May 8, 2019

# Compare log2(Sense) and log2(Antisense) for each gene, mutant v WT

# Use: ./antisense_geneLength.R mut

# Load libraries
library(ggplot2)
library(svglite)

# Get arguments
args <- commandArgs(trailingOnly = TRUE)
mut <- args[1]

# Read table
pos_dat <- read.table(paste("./CodingAntisense/", mut, ".geneNETseqRatio.pos.bed", sep=""), stringsAsFactors = FALSE)
neg_dat <- read.table(paste("./CodingAntisense/", mut, ".geneNETseqRatio.neg.bed", sep=""), stringsAsFactors = FALSE)
mut_dat <- rbind(pos_dat, neg_dat)
mut_dat$Mutant = mut

wt_pos_dat = read.table("./CodingAntisense/wt.geneNETseqRatio.pos.bed", stringsAsFactors = FALSE)
wt_neg_dat = read.table("./CodingAntisense/wt.geneNETseqRatio.neg.bed", stringsAsFactors = FALSE)
wt_dat <- rbind(wt_pos_dat, wt_neg_dat)
wt_dat$Mutant = "wt"

data = rbind(mut_dat, wt_dat)
colnames(data) <- c("Chr", "Start", "End", "Gene", "Value", "Strand", "Sense", "Antisense", "Ratio", "Mutant")

# Make bins for gene length
data$GeneLength = data$End - data$Start

# Correlation
c = cor(data[data$Mut == mut, 'Antisense'], data[data$Mut == mut, 'GeneLength'], use = "pairwise.complete.obs")
cat(paste("There is a Pearson correlation of ", round(c, 2), " between antisense RPKM and gene length\n", sep=""))

# Perform correlation test
cor.test(data[data$Mut == mut, 'Antisense'], data[data$Mut == mut, 'GeneLength'])

# Bin
data$LengthBin = "zzz"
data[data$GeneLength <= 500, 'LengthBin'] = "a"
data[data$GeneLength > 500 & data$GeneLength <= 1000, 'LengthBin'] = "b"
data[data$GeneLength > 1000 & data$GeneLength <= 2000, 'LengthBin'] = "c"
data[data$GeneLength > 2000 & data$GeneLength <= 3000, 'LengthBin'] = "d"
data[data$GeneLength > 3000, 'LengthBin'] = "e"

# Create data frame for plotting
df = data.frame("Bin" = rep(c("a", "b", "c", "d", "e"), 2), Mean = rep(0, 10), SE = rep(0, 10), Mut = rep(c("Mut", "Wt"), each = 5), stringsAsFactors = FALSE)
df[df$Bin == "a" & df$Mut == "Mut", 'Mean'] = mean(data[data$LengthBin == "a" & data$Mut != "wt", 'Antisense'], na.rm = T)
df[df$Bin == "a" & df$Mut == "Mut", 'SE'] = sd(data[data$LengthBin == "a" & data$Mut != "wt", 'Antisense'], na.rm = T) / sqrt(nrow(data[data$LengthBin == "a" & data$Mut != "wt",]))
df[df$Bin == "b" & df$Mut == "Mut", 'Mean'] = mean(data[data$LengthBin == "b" & data$Mut != "wt", 'Antisense'], na.rm = T)
df[df$Bin == "b" & df$Mut == "Mut", 'SE'] = sd(data[data$LengthBin == "b" & data$Mut != "wt", 'Antisense'], na.rm = T) / sqrt(nrow(data[data$LengthBin == "b" & data$Mut != "wt",]))
df[df$Bin == "c" & df$Mut == "Mut", 'Mean'] = mean(data[data$LengthBin == "c" & data$Mut != "wt", 'Antisense'], na.rm = T)
df[df$Bin == "c" & df$Mut == "Mut", 'SE'] = sd(data[data$LengthBin == "c" & data$Mut != "wt", 'Antisense'], na.rm = T) / sqrt(nrow(data[data$LengthBin == "c" & data$Mut != "wt",]))
df[df$Bin == "d" & df$Mut == "Mut", 'Mean'] = mean(data[data$LengthBin == "d" & data$Mut != "wt", 'Antisense'], na.rm = T)
df[df$Bin == "d" & df$Mut == "Mut", 'SE'] = sd(data[data$LengthBin == "d" & data$Mut != "wt", 'Antisense'], na.rm = T) / sqrt(nrow(data[data$LengthBin == "d" & data$Mut != "wt",]))
df[df$Bin == "e" & df$Mut == "Mut", 'Mean'] = mean(data[data$LengthBin == "e" & data$Mut != "wt", 'Antisense'], na.rm = T)
df[df$Bin == "e" & df$Mut == "Mut", 'SE'] = sd(data[data$LengthBin == "e" & data$Mut != "wt", 'Antisense'], na.rm = T) / sqrt(nrow(data[data$LengthBin == "e" & data$Mut != "wt",]))

df[df$Bin == "a" & df$Mut == "Wt", 'Mean'] = mean(data[data$LengthBin == "a" & data$Mut == "wt", 'Antisense'], na.rm = T)
df[df$Bin == "a" & df$Mut == "Wt", 'SE'] = sd(data[data$LengthBin == "a" & data$Mut == "wt", 'Antisense'], na.rm = T) / sqrt(nrow(data[data$LengthBin == "a" & data$Mut == "wt",]))
df[df$Bin == "b" & df$Mut == "Wt", 'Mean'] = mean(data[data$LengthBin == "b" & data$Mut == "wt", 'Antisense'], na.rm = T)
df[df$Bin == "b" & df$Mut == "Wt", 'SE'] = sd(data[data$LengthBin == "b" & data$Mut == "wt", 'Antisense'], na.rm = T) / sqrt(nrow(data[data$LengthBin == "b" & data$Mut == "wt",]))
df[df$Bin == "c" & df$Mut == "Wt", 'Mean'] = mean(data[data$LengthBin == "c" & data$Mut == "wt", 'Antisense'], na.rm = T)
df[df$Bin == "c" & df$Mut == "Wt", 'SE'] = sd(data[data$LengthBin == "c" & data$Mut == "wt", 'Antisense'], na.rm = T) / sqrt(nrow(data[data$LengthBin == "c" & data$Mut == "wt",]))
df[df$Bin == "d" & df$Mut == "Wt", 'Mean'] = mean(data[data$LengthBin == "d" & data$Mut == "wt", 'Antisense'], na.rm = T)
df[df$Bin == "d" & df$Mut == "Wt", 'SE'] = sd(data[data$LengthBin == "d" & data$Mut == "wt", 'Antisense'], na.rm = T) / sqrt(nrow(data[data$LengthBin == "d" & data$Mut == "wt",]))
df[df$Bin == "e" & df$Mut == "Wt", 'Mean'] = mean(data[data$LengthBin == "e" & data$Mut == "wt", 'Antisense'], na.rm = T)
df[df$Bin == "e" & df$Mut == "Wt", 'SE'] = sd(data[data$LengthBin == "e" & data$Mut == "wt", 'Antisense'], na.rm = T) / sqrt(nrow(data[data$LengthBin == "e" & data$Mut == "wt",]))


# Graph
p <- ggplot(df, aes(x=Bin, y=Mean, group = Mut, color = Mut)) + 
  geom_line() + 
  geom_point() + 
  geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=.2, position=position_dodge(0.02)) +
  scale_colour_manual(values=c("Mut" = "black", "Wt" = "grey50")) + 
  theme_bw() + 
  theme(legend.position = "none") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(plot = p, filename = paste(mut, "vWt_antisenseLength.svg", sep=""), height = 5, width = 5)
