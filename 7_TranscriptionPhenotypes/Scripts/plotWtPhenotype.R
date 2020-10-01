#!/n/app/R/3.6.1/bin/Rscript

# Plot number of strains gene is DE in v. transcription phenotypes

# Written by: K. Lachance
# Date: August 28, 2020

# Use: ./plotWtPhenotype.R

library(ggplot2)
library(svglite)

# Get data
DEgenesUP = read.table("DEgenesUP.allStrainsCount.txt", stringsAsFactors = FALSE)
colnames(DEgenesUP) <- c("Count", "Gene")

DEgenesDOWN = read.table("DEgenesDOWN.allStrainsCount.txt", stringsAsFactors = FALSE)
colnames(DEgenesDOWN) <- c("Count", "Gene")

ASdata = read.table("../3_Antisense/CodingAntisense/ALL.geneNETseqRatio.bed", stringsAsFactors = FALSE)
colnames(ASdata) <- c("Chr", "Start", "End", "Gene", "Skip", "Strand", "Sense", "Anti", "Ratio", "Mut")
ASdata = ASdata[ASdata$Mut == "wt",]

TSSdata = read.table("../4_PausingIndex/wt_TSS.PI.bed", stringsAsFactors = FALSE)
colnames(TSSdata) = c("Chr", "Start1", "End1", "Start2", "End2", "Gene", "Skip", "Strand", "ROIcount", "notROIcount", "ROInorm", "notROInorm", "PI", "Region", "Mut")

# Find number of DE genes in WT
genes = unique(c(ASdata$Gene, TSSdata$Gene))
nGenes = length(genes)

# Make dataframe
dfUP = data.frame(Gene = genes, Count = 0, ASratio = 0, TSSpi = 0, PApi = 0, ANTIpi = 0, Dir = "UP", stringsAsFactors = FALSE)
dfDOWN = data.frame(Gene = genes, Count = 0, ASratio = 0, TSSpi = 0, PApi = 0, ANTIpi = 0, Dir = "DOWN", stringsAsFactors = FALSE)

# Loop through genes
cat("Working on up-regulated genes...\n")
for (i in 1:nGenes) {

    if (i %% 1000 == 0) { cat(paste("Done with ", i, " / ", nGenes, " genes!\n", sep="")) }
    
    g = dfUP$Gene[i]
    if (g %in% DEgenesUP$Gene) { dfUP$Count[i] = DEgenesUP[DEgenesUP$Gene == g, 'Count'] }

    if (g %in% ASdata$Gene) { dfUP$ASratio[i] = ASdata[ASdata$Gene == g, 'Ratio'] }
    if (g %in% TSSdata$Gene) { dfUP$TSSpi[i] = TSSdata[TSSdata$Gene == g, 'PI'] }

}

cat("Working on down-regulated genes...\n")
for (i in 1:nGenes) {

    if (i %% 1000 == 0) { cat(paste("Done with ", i, " / ", nGenes, " genes!\n", sep="")) }

    g = dfDOWN$Gene[i]
    if (g %in% DEgenesDOWN$Gene) { dfDOWN$Count[i] = DEgenesDOWN[DEgenesDOWN$Gene == g, 'Count'] }

    if (g %in% ASdata$Gene) { dfDOWN$ASratio[i] = ASdata[ASdata$Gene == g, 'Ratio'] }
    if (g %in% TSSdata$Gene) { dfDOWN$TSSpi[i] = TSSdata[TSSdata$Gene == g, 'PI'] }
}

df = rbind(dfUP, dfDOWN)

df$Cat = NA
df[df$Count == 0 & df$Dir == "UP", 'Cat'] = "UpLess"
df[df$Count == 0 & df$Dir == "DOWN", 'Cat'] = "DownLess"
df[df$Count >= 8 & df$Dir == "UP", 'Cat'] = "Up" # Frequently-regulated genes are those DE in ≥ 8 strains
df[df$Count >= 8 & df$Dir == "DOWN", 'Cat'] = "Down"

DF = df[!is.na(df$Cat),]
df = DF

# T-test

DF_ASratio = df[,c('Cat', 'ASratio')]
DF_ASratio$Log = log2(DF_ASratio$ASratio)
DF_ASratio = DF_ASratio[is.finite(DF_ASratio$Log), ]

DF_TSSpi = df[,c('Cat', 'TSSpi')]
DF_TSSpi$Log = log2(DF_TSSpi$TSSpi)
DF_TSSpi = DF_TSSpi[is.finite(DF_TSSpi$Log), ]

cat("\n AS RATIO \n")
t.test(DF_ASratio[DF_ASratio$Cat == "UpLess", 'Log'], DF_ASratio[DF_ASratio$Cat == "Up", 'Log'])
t.test(DF_ASratio[DF_ASratio$Cat == "DownLess", 'Log'], DF_ASratio[DF_ASratio$Cat == "Down", 'Log'])

cat("\n TSS PI \n")
t.test(DF_TSSpi[DF_TSSpi$Cat == "UpLess", 'Log'], DF_TSSpi[DF_TSSpi$Cat == "Up", 'Log'])
t.test(DF_TSSpi[DF_TSSpi$Cat == "DownLess", 'Log'], DF_TSSpi[DF_TSSpi$Cat == "Down", 'Log'])

# Plot

df$Cat = factor(df$Cat, levels = c("UpLess", "Up", "DownLess", "Down"))

p1 <- ggplot(df, aes(x = Cat, y = log2(ASratio), group = Cat, color = Cat)) + 
   geom_boxplot(outlier.shape = NA) + 
   scale_color_manual(values = c("UpLess" = "#87A6C4", "Up" = "#1874cd", "DownLess" = "#C49E96", "Down" = "#cd5b45")) +
   scale_x_discrete(labels = c("Up (=0)", "Up (>8)", "Down (=0)", "Down (>8)")) +
   theme_bw() + 
   xlab("# Strains Gene is DE") + ylab("log2(AS : S Ratio)") + 
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
   theme(legend.position = "none")

ggsave(plot = p1, filename = "DEgenes_inWT_AntisenseSenseRatio.svg", height = 3, width = 10)

p2 <- ggplot(df, aes(x = Cat, y = log2(TSSpi), group = Cat, color = Cat)) +
   geom_boxplot(outlier.shape = NA) +	     
   scale_color_manual(values = c("UpLess" = "#87A6C4", "Up" = "#1874cd", "DownLess" = "#C49E96", "Down" = "#cd5b45")) +
   scale_x_discrete(labels = c("Up (=0)", "Up (>10)", "Down (=0)", "Down (>10)")) +
   theme_bw() +
   xlab("# Strains Gene is DE") + ylab("log2(TSS PI)") +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
   theme(legend.position = "none")

ggsave(plot = p2, filename = "DEgenes_inWT_TSSpi.svg", height = 3, width = 10)

