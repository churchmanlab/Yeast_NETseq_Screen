#!/n/app/R/4.0.1/bin/Rscript

# Show overlap of pauses with gene regions (start, mid, end) compared to control

# Use: ./Scripts/plotOverlap.R


# Load libraries
library(ggplot2)
library(reshape2)
library(svglite)
library(DescTools)

# Get arguments
args <- commandArgs(trailingOnly = TRUE)

# Load data
df = read.table("./geneOverlap.txt", stringsAsFactors = FALSE)
colnames(df) = c("Mut", "Start", "Mid", "End", "Total")

# Loop through each mut
DF = data.frame(matrix(nrow = nrow(df), ncol = 9, 0), stringsAsFactors = FALSE)
colnames(DF) = c("pStart", "pMid", "pEnd", "lowStart", "lowMid", "lowEnd", "highStart", "highMid", "highEnd")

# Calculate proportion and CI
for (i in 1:nrow(df)) {

    observed = c(df$Start[i], df$Mid[i], df$End[i])
    CI = MultinomCI(observed, conf.level = 0.95, method = "goodman")	#Methods: "sisonglaz", "cplus1", "goodman"
    DF[i,] = CI
}

pDF = DF[,1:3]
pDF$Mut = df$Mut

# Melt, add color categories
M = melt(pDF)
colnames(M) = c("Mut", "Pos", "Value")

M$Cat = as.character(M$Pos)
M[M$Mut == "zzzCNTL" & M$Pos == "pStart", 'Cat'] = "cStart"
M[M$Mut	== "zzzCNTL" & M$Pos == "pMid", 'Cat'] = "cMid"
M[M$Mut	== "zzzCNTL" & M$Pos == "pEnd", 'Cat'] = "cEnd"

# Reorder
ord = pDF$Mut[order(pDF$pStart)] # pStart pEnd
M$Mut <- factor(M$Mut, levels = ord)
M$Mut <- relevel(M$Mut, "zzzCNTL")
M$Mut <- factor(M$Mut, levels=rev(levels(M$Mut)))

M$Cat <- factor(M$Cat, levels = c("pEnd", "pMid", "pStart", "cEnd", "cMid", "cStart"))

# Now deal with the error bars
eDF = DF
eDF$Mut = df$Mut

# Add so that bars appear in correct place
eDF$correct_lowMid = eDF$pStart + eDF$lowMid
eDF$correct_highMid = eDF$pStart + eDF$highMid
eDF$correct_lowEnd = eDF$pStart + eDF$pMid + eDF$lowEnd
eDF$correct_highEnd = eDF$pStart  + eDF$pMid + eDF$highEnd

# Collect corrected errors
E = data.frame(Mut = rep(eDF$Mut, 3), ErrorLow = c(eDF$correct_lowEnd, eDF$correct_lowMid, eDF$lowStart), ErrorHigh = c(eDF$correct_highEnd, eDF$correct_highMid, eDF$highStart), stringsAsFactors = FALSE)

# Reorder
E$Mut <- factor(E$Mut, levels=rev(levels(M$Mut)))

# Plot
p <- ggplot(M, aes(x = Mut, y = Value, fill = Cat)) +
  geom_bar(position = position_stack(), stat = "identity", width = .8) +
  scale_fill_manual(values = c("pStart" = "#356c5a", "pMid" = "#76eec6", "pEnd" = "#ccffee", "cStart" = "grey25", "cMid" = "grey50", "cEnd" = "grey75")) +
  geom_errorbar(aes(ymin=E$ErrorLow, ymax=E$ErrorHigh), width=.2, position="identity") + 
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none", panel.grid = element_blank(),panel.background = element_blank())

ggsave(filename="PauseGeneOverlap.pdf",device="pdf", plot=p, width=6, height=6)




