#!/n/app/R/3.6.1/bin/Rscript

# Make vector for each deletion strain for downstream PCA / clustering

# Written by: K. Lachance
# Date: July 14, 2019

# Use: ./makeMutVector.R AllPauseLoci.bed mut

# Get input parameters
args <- commandArgs(trailingOnly = TRUE)
allPause_file <- args[1]
mut <- args[2]

allPause <- read.table(allPause_file, stringsAsFactors = FALSE)
colnames(allPause) <- c("Chr", "Start", "End", "Strand")
mutPause <- read.table(paste(mut, "_PASS_pauseScore.txt", sep=""), stringsAsFactors = FALSE)
colnames(mutPause) <- c("Chr", "Start", "End", "Score1", "Score2", "Strand", "Norm1", "Norm2")
posPause <- read.table(paste(mut, ".possiblePauses.bed", sep=""), stringsAsFactors = FALSE)
colnames(posPause) <- c("Chr", "Start", "End", "Strand")

# Make new data frame to hold final vector
out <- allPause
out$Vector <- 0

# Loop through each position in posPause and record pause in output
for (i in 1:nrow(posPause)) {
    out[out$Chr == posPause$Chr[i] & out$Start == posPause$Start[i] & out$End == posPause$End[i] & out$Strand == posPause$Strand[i], 'Vector'] = -1
}


# Loop through each position in mutPause and record pause in output
for (i in 1:nrow(mutPause)) {
    out[out$Chr == mutPause$Chr[i] & out$Start == mutPause$Start[i] & out$End == mutPause$End[i] & out$Strand == mutPause$Strand[i], 'Vector'] = 1
}

v = c(mut, out$V)

write.table(t(v), file="mutPauseVectors.txt", append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
