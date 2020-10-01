#!/n/app/R/3.6.1/bin/Rscript

# Plot boxplot of phenotypes, split by regulation status, split by wildtype value

# Written by: K. Lachance
# Date: August 28, 2020

# Use: ./plotPhenotypeSplitByWildtype.R UP/DOWN

library(ggplot2)
library(svglite)

args <- commandArgs(trailingOnly = TRUE)
DIR <- args[1]

# Get data
ASSdata = read.table(paste("DE_AntisenseRatio.", DIR, ".txt", sep=""), stringsAsFactors = FALSE)
colnames(ASSdata) = c("Gene", "WTvalue", "MUTvalue", "Mut", "Cat")

TSSdata = read.table(paste("DE_TssPI.", DIR, ".txt", sep=""), stringsAsFactors = FALSE)
colnames(TSSdata) = c("Gene", "WTvalue", "MUTvalue", "Mut", "Cat")

# Create bins
ASSdata$Bin = "A"
ASSdata[log2(ASSdata$WTvalue) >= -4 & log2(ASSdata$WTvalue) < -2, 'Bin'] = "B"
ASSdata[log2(ASSdata$WTvalue) >= -2 & log2(ASSdata$WTvalue) < 1, 'Bin'] = "C"
ASSdata[log2(ASSdata$WTvalue) >= 1, 'Bin'] = "D"

table(ASSdata[,c('Bin', 'Cat')])

TSSdata$Bin = "A"
TSSdata[TSSdata$WTvalue >= 0.5 & TSSdata$WTvalue < 1.5, 'Bin'] = "B"
TSSdata[TSSdata$WTvalue	>= 1.5 & TSSdata$WTvalue < 2.5, 'Bin'] = "C"
TSSdata[TSSdata$WTvalue	>= 2.5, 'Bin'] = "D"

table(TSSdata[,c('Bin', 'Cat')])

# T test

DF_ASratio = ASSdata[,c('Bin', 'Cat', 'WTvalue')]
DF_ASratio$Log = log2(DF_ASratio$WTvalue)
DF_ASratio = DF_ASratio[is.finite(DF_ASratio$Log), ]

DF_TSSpi = TSSdata
DF_TSSpi$Log = TSSdata$WTvalue

cat("\n AS RATIO \n")
t.test(DF_ASratio[DF_ASratio$Bin == "A" & DF_ASratio$Cat == "NO", 'Log'], DF_ASratio[DF_ASratio$Bin == "A" & DF_ASratio$Cat == "YES", 'Log'])
t.test(DF_ASratio[DF_ASratio$Bin == "B" & DF_ASratio$Cat == "NO", 'Log'], DF_ASratio[DF_ASratio$Bin == "B" & DF_ASratio$Cat == "YES", 'Log'])
t.test(DF_ASratio[DF_ASratio$Bin == "C" & DF_ASratio$Cat == "NO", 'Log'], DF_ASratio[DF_ASratio$Bin == "C" & DF_ASratio$Cat == "YES", 'Log'])
t.test(DF_ASratio[DF_ASratio$Bin == "D" & DF_ASratio$Cat == "NO", 'Log'], DF_ASratio[DF_ASratio$Bin == "D" & DF_ASratio$Cat == "YES", 'Log'])

cat("\n TSS PI \n")
t.test(DF_TSSpi[DF_TSSpi$Bin == "A" & DF_TSSpi$Cat == "NO", 'Log'], DF_TSSpi[DF_TSSpi$Bin == "A" & DF_TSSpi$Cat == "YES", 'Log'])
t.test(DF_TSSpi[DF_TSSpi$Bin ==	"B" & DF_TSSpi$Cat == "NO", 'Log'], DF_TSSpi[DF_TSSpi$Bin == "B" & DF_TSSpi$Cat == "YES", 'Log'])
t.test(DF_TSSpi[DF_TSSpi$Bin ==	"C" & DF_TSSpi$Cat == "NO", 'Log'], DF_TSSpi[DF_TSSpi$Bin == "C" & DF_TSSpi$Cat == "YES", 'Log'])
#t.test(DF_TSSpi[DF_TSSpi$Bin ==	"D" & DF_TSSpi$Cat == "NO", 'Log'], DF_TSSpi[DF_TSSpi$Bin == "D" & DF_TSSpi$Cat == "YES", 'Log'])

# Set colors
if (DIR == "UP") {
   colYES = "#1874cd"
   colNO = "#87A6C4"
} else {
   colYES = "#cd5b45"
   colNO = "#C49E96"
}

# Plot
p1 <- ggplot(ASSdata, aes(x = Bin, y = log2(MUTvalue), color = Cat)) + 
   geom_boxplot(outlier.shape = NA) + 
   scale_color_manual(values = c("YES" = colYES, "NO" = colNO)) +
   coord_cartesian(ylim = c(-10, 6)) + 
   theme_bw() + 
   xlab("Value in Wild-type") + ylab("Value in Deletion Strains") +
   scale_x_discrete(labels = c("A" = "< -5", "B" = "-5 to -3", "C" = "-3 to -1", "D" = "-1 to 1", "E" = "> 1")) +  
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
   theme(legend.position = "none") + 
   ggtitle("AS : S Ratio") + theme(plot.title = element_text(hjust = 0.5)) 

ggsave(plot = p1, filename = paste("DEboxplot_AntisenseRatio.", DIR, ".svg", sep=""), height = 6, width = 10)

p2 <- ggplot(TSSdata, aes(x = Bin, y = MUTvalue, color = Cat)) +
   geom_boxplot(outlier.shape = NA) +
   scale_color_manual(values = c("YES" = colYES, "NO" = colNO)) +
   coord_cartesian(ylim = c(0, 9)) + 
   theme_bw() +
   xlab("Value in Wild-type") + ylab("Value in Deletion Strains") +
   scale_x_discrete(labels = c("A" = "< 0.5", "B" = "0.5 to 1", "C" = "1 to 1.5", "D" = "1.5 to 2", "E" = "2 to 2.5", "F" = "2.5 to 3", "G" = "> 3")) + 
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
   theme(legend.position = "none") +
   ggtitle("TSS PI") + theme(plot.title = element_text(hjust = 0.5))

ggsave(plot = p2, filename = paste("DEboxplot_TSSpi.", DIR, ".svg", sep=""), height = 6, width = 10)
