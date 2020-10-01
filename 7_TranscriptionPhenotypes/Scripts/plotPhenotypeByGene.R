#!/n/app/R/3.6.1/bin/Rscript

# Plot the phenotype of each frequently regulated gene when differentially expressed compared to when not DE

# Written by: K. Lachance
# Date: August 28, 2020

# Use: ./plotPhenotypeByGene.R

library(ggplot2)
library(svglite)
library(dplyr)
library(reshape2)

# Get data
ASSdataUP = read.table("DE_AntisenseRatio.UP.txt", stringsAsFactors = FALSE)
colnames(ASSdataUP) = c("Gene", "WTvalue", "MUTvalue", "Mut", "Cat")
ASSdataUP$Dir = "UP"
ASSdataDOWN = read.table("DE_AntisenseRatio.DOWN.txt", stringsAsFactors = FALSE)
colnames(ASSdataDOWN) = c("Gene", "WTvalue", "MUTvalue", "Mut", "Cat")
ASSdataDOWN$Dir = "DOWN"

TSSdataUP = read.table("DE_TssPI.UP.txt", stringsAsFactors = FALSE)
colnames(TSSdataUP) = c("Gene", "WTvalue", "MUTvalue", "Mut", "Cat")
TSSdataUP$Dir = "UP"
TSSdataDOWN = read.table("DE_TssPI.DOWN.txt", stringsAsFactors = FALSE)
colnames(TSSdataDOWN) = c("Gene", "WTvalue", "MUTvalue", "Mut", "Cat")
TSSdataDOWN$Dir = "DOWN"

# Calculate median values
myMed = function(X) {
      x = X[is.finite(X)]
      return(median(x, na.rm = TRUE))
}

ASS_UP = aggregate(log2(ASSdataUP$MUTvalue), by = list(Gene = ASSdataUP$Gene, Category = ASSdataUP$Cat, Direction = ASSdataUP$Dir), FUN = myMed)
ASS_DOWN = aggregate(log2(ASSdataDOWN$MUTvalue), by = list(Gene = ASSdataDOWN$Gene, Category = ASSdataDOWN$Cat, Direction = ASSdataDOWN$Dir), FUN = myMed)

TSS_UP = aggregate(TSSdataUP$MUTvalue, by = list(Gene = TSSdataUP$Gene, Category = TSSdataUP$Cat, Direction = TSSdataUP$Dir), FUN = myMed)
TSS_DOWN = aggregate(TSSdataDOWN$MUTvalue, by = list(Gene = TSSdataDOWN$Gene, Category = TSSdataDOWN$Cat, Direction = TSSdataDOWN$Dir), FUN = myMed)

# Calculate difference
ASS_med_UP = aggregate(log2(ASSdataUP$MUTvalue), by = list(Gene = ASSdataUP$Gene, Category = ASSdataUP$Cat), FUN = myMed)
dASS_UP = dcast(ASS_med_UP, Gene ~ Category)
dASS_UP$Diff = dASS_UP$YES - dASS_UP$NO
ASSneg_UP = dASS_UP[dASS_UP$Diff < 0, 'Gene']

ASS_med_DOWN = aggregate(log2(ASSdataDOWN$MUTvalue), by = list(Gene = ASSdataDOWN$Gene, Category = ASSdataDOWN$Cat), FUN = myMed)
dASS_DOWN = dcast(ASS_med_DOWN, Gene ~ Category)
dASS_DOWN$Diff = dASS_DOWN$YES - dASS_DOWN$NO
ASSneg_DOWN = dASS_DOWN[dASS_DOWN$Diff < 0, 'Gene']

TSS_med_UP = aggregate(TSSdataUP$MUTvalue, by = list(Gene = TSSdataUP$Gene, Category = TSSdataUP$Cat), FUN = myMed)
dTSS_UP = dcast(TSS_med_UP, Gene ~ Category)
dTSS_UP$Diff = dTSS_UP$YES - dTSS_UP$NO
TSSneg_UP = dTSS_UP[dTSS_UP$Diff < 0, 'Gene']

TSS_med_DOWN = aggregate(TSSdataDOWN$MUTvalue, by = list(Gene = TSSdataDOWN$Gene, Category = TSSdataDOWN$Cat), FUN = myMed)
dTSS_DOWN = dcast(TSS_med_DOWN, Gene ~ Category)
dTSS_DOWN$Diff = dTSS_DOWN$YES - dTSS_DOWN$NO
TSSneg_DOWN = dTSS_DOWN[dTSS_DOWN$Diff < 0, 'Gene']

# Add color categories
ASS_UP$COLOR = "UpYes"
ASS_UP[ASS_UP$Cat == "NO", 'COLOR'] = "UpNo"

ASS_DOWN$COLOR = "DownYes"
ASS_DOWN[ASS_DOWN$Cat == "NO", 'COLOR'] = "DownNo"

TSS_UP$COLOR = "UpYes"
TSS_UP[TSS_UP$Cat == "NO", 'COLOR'] = "UpNo"

TSS_DOWN$COLOR = "DownYes"
TSS_DOWN[TSS_DOWN$Cat == "NO", 'COLOR'] = "DownNo"

# Relevel for plotting
ASS_UP$Gene = factor(ASS_UP$Gene, levels = dASS_UP$Gene[order(dASS_UP$Diff)])
ASS_DOWN$Gene = factor(ASS_DOWN$Gene, levels = dASS_DOWN$Gene[order(dASS_DOWN$Diff)])

TSS_UP$Gene = factor(TSS_UP$Gene, levels = dTSS_UP$Gene[order(dTSS_UP$Diff)])
TSS_DOWN$Gene = factor(TSS_DOWN$Gene, levels = dTSS_DOWN$Gene[order(dTSS_DOWN$Diff)])

# Set colors
UP_YES = "#1874cd"
UP_NO = "#87A6C4"
DOWN_YES = "#cd5b45"
DOWN_NO = "#C49E96"

# Plot
p1 <- ggplot(ASS_UP, aes(x = Gene, y = x, color = COLOR, group = Gene)) + 
   geom_line() + 
   geom_point() + 
   scale_color_manual(values = c("UpYes" = UP_YES, "UpNo" = UP_NO, "DownYes" = DOWN_YES, "DownNo" = DOWN_NO)) +
   theme_bw() + 
#   coord_cartesian(ylim = c(-8, 8)) + 
   xlab("Gene") + ylab("Value in Deletion Strains") +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
   theme(legend.position = "none") + 
   ggtitle("AS : S Ratio") + theme(plot.title = element_text(hjust = 0.5)) 

ggsave(plot = p1, filename = "DEbyGeneMed_AntisenseRatio.UP.svg", height = 3, width = 12)

p2 <- ggplot(ASS_DOWN, aes(x = Gene, y = x, color = COLOR, group = Gene)) +
   geom_line() +
   geom_point() +
   scale_color_manual(values = c("UpYes" = UP_YES, "UpNo" = UP_NO, "DownYes" = DOWN_YES, "DownNo" = DOWN_NO)) +
   theme_bw() +
#   coord_cartesian(ylim = c(-8, 8)) +
   xlab("Gene") + ylab("Value in Deletion Strains") +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
   theme(legend.position = "none") +
   ggtitle("AS : S Ratio") + theme(plot.title = element_text(hjust = 0.5)) 

ggsave(plot = p2, filename = "DEbyGeneMed_AntisenseRatio.DOWN.svg", height = 3, width = 12)

p3 <- ggplot(TSS_UP, aes(x = Gene, y = x, color = COLOR, group = Gene)) +
   geom_line() +
   geom_point() +
   scale_color_manual(values = c("UpYes" = UP_YES, "UpNo" = UP_NO, "DownYes" = DOWN_YES, "DownNo" = DOWN_NO)) +
   theme_bw() +
#   coord_cartesian(ylim	= c(0, 15)) +	
   xlab("Gene") + ylab("Value in Deletion Strains") +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
   theme(legend.position = "none") +
   ggtitle("TSS PI") + theme(plot.title = element_text(hjust = 0.5)) 

ggsave(plot = p3, filename = "DEbyGeneMed_TSSpi.UP.svg", height = 3, width = 12)

p4 <- ggplot(TSS_DOWN, aes(x = Gene, y = x, color = COLOR, group = Gene)) +
   geom_line() +
   geom_point() +
   scale_color_manual(values = c("UpYes" = UP_YES, "UpNo" = UP_NO, "DownYes" = DOWN_YES, "DownNo" = DOWN_NO)) +
   theme_bw() +
#   coord_cartesian(ylim = c(0, 15)) +
   xlab("Gene") + ylab("Value in Deletion Strains") +
   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
   theme(legend.position = "none") +
   ggtitle("TSS PI") + theme(plot.title = element_text(hjust = 0.5)) 

ggsave(plot = p4, filename = "DEbyGeneMed_TSSpi.DOWN.svg", height = 3, width = 12)