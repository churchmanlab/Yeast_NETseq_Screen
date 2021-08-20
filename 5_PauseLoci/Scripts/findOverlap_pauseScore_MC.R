#!/n/app/R/4.0.1/bin/Rscript

# Find overlapping pauses while maintaining scores (for IDR)

# Written by: K. Lachance
# Date: Feb 5, 2020
# Update by: M. Couvillion
# Date: 5/2021

# Use: ./findOverlap_pauseScore_MC.R mut w covT

# Get input parameters
args <- commandArgs(trailingOnly = TRUE)
MUT <- args[1]
wind <- args[2]
covT <- args[3]

# Load data
pos1 = read.table(paste(MUT, "-1_cov",covT, ".pausesNB_T3.w",wind,".pos.bed", sep=""), stringsAsFactors = FALSE)
colnames(pos1) = c("Chr", "Start", "End", "Score1")
pos2 = read.table(paste(MUT, "-2_cov",covT, ".pausesNB_T3.w",wind,".pos.bed", sep=""), stringsAsFactors = FALSE)
colnames(pos2) = c("Chr", "Start", "End", "Score2")

neg1 = read.table(paste(MUT, "-1_cov",covT, ".pausesNB_T3.w",wind,".neg.bed", sep=""), stringsAsFactors = FALSE)
colnames(neg1) = c("Chr", "Start", "End", "Score1")
neg2 = read.table(paste(MUT, "-2_cov",covT, ".pausesNB_T3.w",wind,".neg.bed", sep=""), stringsAsFactors = FALSE)
colnames(neg2) = c("Chr", "Start", "End", "Score2")

# Combine pauses by strand
POS = unique(merge(pos1, pos2, by=c("Chr", "Start", "End")))
NEG = unique(merge(neg1, neg2, by=c("Chr", "Start", "End")))

# Print information about overlap
n1 = nrow(pos1) + nrow(neg1)
n2 = nrow(pos2) + nrow(neg2)
nComb = nrow(POS) + nrow(NEG)

tmp = data.frame(MUT, n1, n2, nComb, round((nComb / (n1 + n2 - nComb)) * 100, 2), stringsAsFactors = FALSE)
write.table(tmp, paste0("IDR_reproducibility_cov", covT, "_tmp.txt"), append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Add strand information
POS$Strand = "+"
NEG$Strand = "-"

# Combine strands
PAUSES = rbind(POS, NEG)

# Normalize by RPKM
readCount = read.table("./readTotal.txt", stringsAsFactors = FALSE)
colnames(readCount) = c("Mut", "Count")

COUNT1 = readCount[readCount$Mut == paste0(MUT, "-1"), 'Count']
COUNT2 = readCount[readCount$Mut == paste0(MUT, "-2"), 'Count']

PAUSES$NormScore1 = (PAUSES$Score1 / COUNT1) * 1000000
PAUSES$NormScore2 = (PAUSES$Score2 / COUNT2) * 1000000

# Write output table
write.table(PAUSES, paste0(MUT, "_OverlapScore_cov",covT,".txt"), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

