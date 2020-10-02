#!/n/app/R/3.6.1/bin/Rscript

# Consolidate overlap and combine into single file

# Written by: K. Lachance
# Date: May 10, 2020

# Use: ./calcPI.R mut region

# Load libraries
suppressMessages(library(dplyr))

# Get arguments
args <- commandArgs(trailingOnly = TRUE)
mut <- args[1]
label <- args[2]

# Read in tables
ROI_pos <- read.table(paste(mut, ".tmp.pos.bed", sep=""), stringsAsFactors = FALSE)
colnames(ROI_pos) = c("GeneChr", "GeneStart", "GeneEnd", "GeneName", "Skip", "GeneStrand", "NETChr", "NETStart", "NETEnd", "NETValue")
ROI_neg <- read.table(paste(mut, ".tmp.neg.bed", sep=""), stringsAsFactors = FALSE)
colnames(ROI_neg) = c("GeneChr", "GeneStart", "GeneEnd", "GeneName", "Skip", "GeneStrand", "NETChr", "NETStart", "NETEnd", "NETValue")

NOT_pos <- read.table(paste(mut, ".tmpNOT.pos.bed", sep=""), stringsAsFactors = FALSE)
colnames(NOT_pos) = c("GeneChr", "GeneStart", "GeneEnd", "GeneName", "Skip", "GeneStrand", "NETChr", "NETStart", "NETEnd", "NETValue")
NOT_neg <- read.table(paste(mut, ".tmpNOT.neg.bed", sep=""), stringsAsFactors = FALSE)
colnames(NOT_neg) = c("GeneChr", "GeneStart", "GeneEnd", "GeneName", "Skip", "GeneStrand", "NETChr", "NETStart", "NETEnd", "NETValue")

# Consolidate data by gene
ROIp = ROI_pos[,1:6] %>% distinct(GeneName, .keep_all = TRUE)
counts = aggregate(ROI_pos$NETValue, by=list(GeneName=ROI_pos$GeneName), FUN=sum)
Rp = merge(x = ROIp, y = counts, by = "GeneName", all.x = TRUE)
RP = Rp[,c(2,3,4,1,5,6,7)]
colnames(RP) = c("Chr", "Start", "End", "Gene", "Skip", "Strand", "ReadCount")

ROIn = ROI_neg[,1:6] %>% distinct(GeneName, .keep_all = TRUE)
counts = aggregate(ROI_neg$NETValue, by=list(GeneName=ROI_neg$GeneName), FUN=sum)
Rn = merge(x = ROIn, y = counts, by = "GeneName", all.x = TRUE)
RN = Rn[,c(2,3,4,1,5,6,7)]
colnames(RN) = c("Chr",	"Start", "End",	"Gene",	"Skip",	"Strand", "ReadCount")

NOTp = NOT_pos[,1:6] %>% distinct(GeneName, .keep_all = TRUE)
counts = aggregate(NOT_pos$NETValue, by=list(GeneName=NOT_pos$GeneName), FUN=sum)
Np = merge(x = NOTp, y = counts, by = "GeneName", all.x = TRUE)
NP = Np[,c(2,3,4,1,5,6,7)]
colnames(NP) = c("Chr",	"Start", "End",	"Gene",	"Skip",	"Strand", "ReadCount")

NOTn = NOT_neg[,1:6] %>% distinct(GeneName, .keep_all = TRUE)
counts = aggregate(NOT_neg$NETValue, by=list(GeneName=NOT_neg$GeneName), FUN=sum)
Nn = merge(x = NOTn, y = counts, by = "GeneName", all.x = TRUE)
NN = Nn[,c(2,3,4,1,5,6,7)]
colnames(NN) = c("Chr",	"Start", "End",	"Gene",	"Skip",	"Strand", "ReadCount")

# Combine
ROI = rbind(RP, RN)
NOT = rbind(NP, NN)

# Merge
tmp = merge(x = ROI, y = NOT, by = "Gene")
out = tmp[,c(2,3,4,9,10,1,5,6,7,13)]
colnames(out) = c("Chr", "ROIStart", "ROIEnd", "NOTStart", "NOTEnd", "Gene", "Skip", "Strand", "ROIReadCount", "NOTReadCount")

# Calculate normalized value
out$ROIReadNorm = out$ROIReadCount / (out$ROIEnd - out$ROIStart)
out$NOTReadNorm = out$NOTReadCount / (out$NOTEnd - out$NOTStart)

# Calculate PI
out$PI = out$ROIReadNorm / out$NOTReadNorm

# Add labels
out$Label = paste(label)
out$Mut = paste(mut)

write.table(out, paste(mut, "_", label, ".PI.bed", sep=""), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")