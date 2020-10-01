#!/n/app/R/3.6.1/bin/Rscript

# Plot number of strains gene is DE in v. transcription phenotypes

# Written by: K. Lachance
# Date: August 28, 2020

# Use: ./plotExampleGene.R gene UP/DOWN

library(ggplot2)
library(svglite)
library(dplyr)
library(reshape2)

# Read in parameters (data)
args <- commandArgs(trailingOnly = TRUE)
GENE <- args[1]
DIR <- args[2]

# Get data
ASSdata = read.table(paste("DE_AntisenseRatio.", DIR, ".txt", sep=""), stringsAsFactors = FALSE)
colnames(ASSdata) = c("Gene", "WTvalue", "MUTvalue", "Mut", "Cat")

TSSdata = read.table(paste("DE_TssPI.", DIR, ".txt", sep=""), stringsAsFactors = FALSE)
colnames(TSSdata) = c("Gene", "WTvalue", "MUTvalue", "Mut", "Cat")

# Get gene
ASS = ASSdata[ASSdata$Gene == GENE, ]
TSS = TSSdata[TSSdata$Gene == GENE, ]

# Calculate medians
mYES_ASS = median(log2(ASS[ASS$Cat == "YES", 'MUTvalue']), na.rm = TRUE)
mNO_ASS = median(log2(ASS[ASS$Cat == "NO", 'MUTvalue']), na.rm = TRUE)
ASS_medians = data.frame(Cat = c("YES", "NO"), Med = c(mYES_ASS, mNO_ASS), stringsAsFactors = FALSE)

mYES_TSS = median(TSS[TSS$Cat == "YES", 'MUTvalue'], na.rm = TRUE)
mNO_TSS = median(TSS[TSS$Cat == "NO", 'MUTvalue'], na.rm = TRUE)
TSS_medians = data.frame(Cat = c("YES", "NO"), Med = c(mYES_TSS, mNO_TSS), stringsAsFactors = FALSE)

# Set colors
if (DIR == "UP") {
   YES = "#1874cd"
   NO = "#87A6C4"
} 
if (DIR == "DOWN") {
   YES = "#cd5b45"
   NO = "#C49E96"
}

# Plot
p1 <- ggplot() + 
   geom_hline(yintercept = log2(ASS$WTvalue[1]), linetype = "dotted", color = "grey50") + 
   geom_jitter(data = ASS, aes(x = Cat, y = log2(MUTvalue), color = Cat), width = 0.2) +
   geom_crossbar(data = ASS_medians, aes(x = Cat, y = Med, color = Cat, ymin = Med, ymax = Med), width = .5) + 
   scale_color_manual(values = c("YES" = YES, "NO" = NO)) + 
   theme_bw() + 
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
   theme(legend.position = "none") 

ggsave(plot = p1, filename = paste("DEboxplotExGene_", GENE, "_AntisenseRatio.", DIR, ".svg", sep=""), height = 3, width = 4)

p2 <- ggplot(TSS, aes(x = Cat, y = MUTvalue, color = Cat)) +
   geom_hline(yintercept = TSS$WTvalue[1], linetype = "dotted", color = "grey50") +
   geom_jitter(data = TSS, aes(x = Cat, y = MUTvalue, color = Cat), width = 0.2) +
   geom_crossbar(data = TSS_medians, aes(x = Cat, y = Med, color = Cat, ymin = Med, ymax = Med), width = .5) +
   scale_color_manual(values = c("YES" = YES, "NO" = NO)) +
   theme_bw() +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
   theme(legend.position = "none") 

ggsave(plot = p2, filename = paste("DEboxplotExGene_", GENE, "_TSSpi.", DIR, ".svg", sep=""), height = 3, width = 4)
