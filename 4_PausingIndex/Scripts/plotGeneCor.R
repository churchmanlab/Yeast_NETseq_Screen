#!/n/app/R/3.6.1/bin/Rscript

# Written by: K. Lachance
# Date: June 9, 2020

# Plot the correlation between PIs for an individual sample

# Use: ./plotGeneCor.R sampleName

library(ggplot2)
library(svglite)

# Read in parameters
args <- commandArgs(trailingOnly = TRUE)
mut <- args[1]

data1 = read.table(paste(mut, "_TSS.PI.bed", sep=""), stringsAsFactors = F)[,c(6,13,14)]
colnames(data1) = c("Gene", "PI", "Region")

data2 = read.table(paste(mut, "_pA.PI.bed", sep=""), stringsAsFactors = F)[,c(6,13,14)]
colnames(data2) = c("Gene", "PI", "Region")

data3 = read.table(paste(mut, "_5pSS.PI.bed", sep=""), stringsAsFactors = F)[,c(6,13,14)]
colnames(data3) = c("Gene", "PI", "Region")

data4 = read.table(paste(mut, "_3pSS.PI.bed", sep=""), stringsAsFactors = F)[,c(6,13,14)]
colnames(data4) = c("Gene", "PI", "Region")

data5 = read.table(paste(mut, "_Antisense.PI.bed", sep=""), stringsAsFactors = F)[,c(6,13,14)]
colnames(data5) = c("Gene", "PI", "Region")

# Merge
df1 = merge(data1, data2, by = "Gene")
df2 = merge(df1, data5, by = "Gene")

df_StartEnd = df2[,c(1,2,4,6)]
colnames(df_StartEnd) = c("Gene", "TSS", "pA", "Anti")

df3 = merge(data4, data5, by = "Gene")

df_Splice = df3[,c(1,2,4)]
colnames(df_Splice) = c("Gene", "Five", "Three")

# Log
df_StartEnd$LogTSS = log2(df_StartEnd$TSS)
df_StartEnd$LogPA = log2(df_StartEnd$pA)
df_StartEnd$LogANTI = log2(df_StartEnd$Anti)

df_StartEnd = df_StartEnd[is.finite(df_StartEnd$LogTSS) & is.finite(df_StartEnd$LogPA) & is.finite(df_StartEnd$LogANTI),]

df_Splice$LogFIVE = log2(df_Splice$Five)
df_Splice$LogTHREE = log2(df_Splice$Three)

df_Splice = df_Splice[is.finite(df_Splice$LogFIVE) & is.finite(df_Splice$LogTHREE),]

# Calculate correlations
c1 = cor(df_StartEnd$LogTSS, df_StartEnd$LogPA, method = "pearson")
cat(paste("There is a Pearson correlation of ", round(c1, 2), " between TSS PI and pA PI.\n", sep=""))

c2 = cor(df_StartEnd$LogTSS, df_StartEnd$LogANTI, method = "pearson")
cat(paste("There is a Pearson correlation of ", round(c2, 2), " between TSS PI and Antisense PI.\n", sep=""))

c3 = cor(df_StartEnd$LogPA, df_StartEnd$LogANTI, method = "pearson")
cat(paste("There is a Pearson correlation of ", round(c3, 2), " between pA PI and Antisense PI.\n", sep=""))

c4 = cor(df_Splice$LogFIVE, df_Splice$LogTHREE, method = "pearson")
cat(paste("There is a Pearson correlation of ", round(c4, 2), " between 5pSS PI and 3pSS PI.\n", sep=""))

# Plot
p1 = ggplot(df_StartEnd, aes(x = LogTSS, y = LogPA)) + 
  geom_point(alpha = 0.4, color = "grey50") + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(filename = paste(mut, "_TSSvPA_PIscatter.svg", sep=""), plot = p1, height = 5, width = 5)

p2 = ggplot(df_StartEnd, aes(x = LogTSS, y = LogANTI)) +
  geom_point(alpha = 0.4, color = "grey50") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(filename = paste(mut, "_TSSvAntisense_PIscatter.svg", sep=""), plot = p2, height = 5, width = 5)

p3 = ggplot(df_StartEnd, aes(x = LogPA, y = LogANTI)) +
  geom_point(alpha = 0.4, color = "grey50") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(filename = paste(mut, "_pAvAntisense_PIscatter.svg", sep=""), plot = p3, height = 5, width = 5)

p4 = ggplot(df_Splice, aes(x = LogFIVE, y = LogTHREE)) +
  geom_point(alpha = 0.4, color = "grey50") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave(filename = paste(mut, "_5pSSv3pSS_PIscatter.svg", sep=""), plot = p4, height = 5, width = 5)