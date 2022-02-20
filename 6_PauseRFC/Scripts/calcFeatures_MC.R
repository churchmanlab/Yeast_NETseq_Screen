#!/n/app/R/4.0.1/bin/Rscript

# Calculate features for each pause using several methods

# Written by: K. Lachance
# Date: April 21, 2020

# Load libraries
library(stringr)
library(TmCalculator)
library(DNAshapeR)

# Get input parameters
args <- commandArgs(trailingOnly = TRUE)
MUT <- args[1]
TYPE <- args[2] # Real / Shuffled

# Command to extract data in each window for each pause
bashCommand='/n/app/bedops/2.4.30/bedextract $1 $2 '

# Read data
if (TYPE == "Real") {
   data = read.table(paste("../5_PauseLoci/", MUT, ".IDRrep.PAUSES.bed", sep=""), stringsAsFactors = FALSE)
   data_seq = read.table(paste("../5_PauseLoci/", MUT, ".IDRrep.pauseWindows.fasta", sep = ""), stringsAsFactors = FALSE)

} else {
   data = read.table(paste("../5_PauseLoci/", MUT, ".IDRrep.SHUFFLE.bed", sep=""), stringsAsFactors = FALSE)
   data_seq = read.table(paste("../5_PauseLoci/", MUT, ".IDRrep.pauseWindowsSHUFFLE.fasta", sep = ""), stringsAsFactors = FALSE)
}
colnames(data) = c("Chr", "Start", "End", "Strand")

data_seq = data_seq[!grepl(">chr", data_seq$V1),]
data$Seq = data_seq

# Function that reverses a string by characters
reverse_chars <- function(string) {
  string_split = strsplit(string, split = "")
  rev_order = nchar(string):1
  reversed_chars = string_split[[1]][rev_order]
  paste(reversed_chars, collapse = "")
} 

# Add features for each pause (proximal nucleotides, DNA melting temp, etc.)
for (i in 1:nrow(data)) {

    if (i %% 100 == 0) {cat(paste("Done with ", i, " / ", nrow(data), " pauses!\n"))}

    # Proximal nucleotides
    if (data$Strand[i] == "+") {
       data$PauseNucA[i] = substr(data$Seq[i], 9, 9)
       data$PauseNucB[i] = substr(data$Seq[i], 10, 10)
       data$PauseNucC[i] = substr(data$Seq[i], 11, 11)
       data$PauseNucD[i] = substr(data$Seq[i], 12, 12)
       data$PauseNucE[i] = substr(data$Seq[i], 13, 13)
       data$PauseNucAB[i] = substr(data$Seq[i], 9, 10)
       data$PauseNucBC[i] = substr(data$Seq[i], 10, 11)
       data$PauseNucCD[i] = substr(data$Seq[i], 11, 12)
       data$PauseNucDE[i] = substr(data$Seq[i], 12, 13)
    } else {
       data$PauseNucA[i] = substr(data$Seq[i], 13, 13)
       data$PauseNucB[i] = substr(data$Seq[i], 12, 12)
       data$PauseNucC[i] = substr(data$Seq[i], 11, 11)
       data$PauseNucD[i] = substr(data$Seq[i], 10, 10)
       data$PauseNucE[i] = substr(data$Seq[i], 9, 9)
       data$PauseNucAB[i] = reverse_chars(substr(data$Seq[i], 12, 13))
       data$PauseNucBC[i] = reverse_chars(substr(data$Seq[i], 11, 12))
       data$PauseNucCD[i] = reverse_chars(substr(data$Seq[i], 10, 11))
       data$PauseNucDE[i] = reverse_chars(substr(data$Seq[i], 9, 10))
    }

    # % CG
    data$pCG[i] = round((str_count(data$Seq[i], "C") + str_count(data$Seq[i], "G")) / 21, 2)

    # Write temporary file for overlaps
    tmp = data.frame(Chr = data$Chr[i], Start = data$Start[i], End = data$End[i], stringsAsFactors = FALSE)
    write.table(tmp, paste(MUT, "_", TYPE, ".tmpCoord.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

    # Overlap with data
#     data_list = c("DMS", "MNase", "h2ak5a", "h3k14ac", "h3k23ac", "h3k36me2", "h3k36me", "h3k4me2", "h3k4me", "h3k79me3", "h3k9ac", "h4k12ac", "h4k20me", "h4k8ac", "h4r3me", "h2as129ph", "h3k18ac", "h3k27ac", "h3k36me3", "h3k4ac", "h3k4me3", "h3k56ac", "h3k79me", "h3s10ph", "h4k16ac", "h4k5ac", "h4r3me2s", "htz1", "Ser2P", "Ser5P", "Ser7P")
    data_list = c("DMS", "MNase", "h2ak5a", "h2as129ph", "h3k14ac", "h3k18ac", "h3k23ac", "h3k27ac", "h3k36me2", "h3k36me3", "h3k36me", "h3k4ac", "h3k4me2", "h3k4me3", "h3k4me", "h3k56ac", "h3k79me3", "h3k79me", "h3k9ac", "h3s10ph", "h4k12ac", "h4k16ac", "h4k20me", "h4k5ac", "h4k8ac", "h4r3me2s", "h4r3me", "htz1", "Ser2P","Ser5P","Ser7P")
    # Loop through each data type
    N = 15 # Number of columns existing before adding these new data types
    for (j in 1:length(data_list)) {
        DATA_TYPE = data_list[j]
	COL = N + j

	temp = read.table(pipe(paste(bashCommand, '../0_Annotations/ChromatinFeatureData/', DATA_TYPE, '.bed ', MUT, '_', TYPE, '.tmpCoord.txt', sep='')), sep="\t", header=FALSE, col.names=c("Chr", "Start", "End", "Value"), colClasses=c('character', rep('numeric', 3)))

    	data[i, COL] = temp$Value[1]
	colnames(data)[COL] <- paste(DATA_TYPE)
    }

    # Nearest nucleosome center
    temp = read.table(pipe(paste('closestBed -t first -a ', MUT, '_', TYPE, '.tmpCoord.txt -b ../0_Annotations/ChromatinFeatureData/nucCenters.bed', sep='')), sep="\t", header=FALSE, col.names=c("PauseChr", "PauseStart", "PauseEnd", "NucChr", "NucStart", "NucEnd"), colClasses=c('character', 'numeric', 'numeric', 'character', 'numeric', 'numeric'))
    data$NucCenter[i] = abs(temp$NucStart[1] - temp$PauseStart[1])

    # Percent into gene body
    # First make new files with transcripts instead of ORFs
    # awk -F "\t" 'BEGIN {OFS="\t"} {$2 = $2 - 200; $3 = $3 + 200; print}' txpts_all_pos.bed > txpts_all_buffer200.pos.bed
    if (data$Strand[i] == "+") { 
       temp = read.table(pipe(paste('intersectBed -wb -a ', MUT, '_', TYPE, '.tmpCoord.txt -b ../0_Annotations/Genes/txpts_all_pos.bed', sep='')), sep="\t", header=FALSE, col.names=c("PauseChr", "PauseStart", "PauseEnd", "GeneChr", "GeneStart", "GeneEnd", "GeneName", "GeneValue", "GeneStrand"), colClasses=c('character', 'numeric', 'numeric', 'character', 'numeric', 'numeric', 'character', 'character', 'character')) 
       data$GenePerc[i] = round(((temp$PauseStart[1] - temp$GeneStart[1]) / (temp$GeneEnd[1] - temp$GeneStart[1])) * 100, 2)
       data$DistFromTSS[i] = temp$PauseStart[1]  - temp$GeneStart[1]
       data$DistFrompA[i] = temp$GeneEnd[1] - temp$PauseStart[1]
    } else {
       temp = read.table(pipe(paste('intersectBed -wb -a ', MUT, '_', TYPE, '.tmpCoord.txt -b ../0_Annotations/Genes/txpts_all_neg.bed', sep='')), sep="\t", header=FALSE, col.names=c("PauseChr", "PauseStart", "PauseEnd", "GeneChr", "GeneStart", "GeneEnd", "GeneName", "GeneValue", "GeneStrand"), colClasses=c('character', 'numeric', 'numeric', 'character', 'numeric', 'numeric', 'character', 'character', 'character'))
       data$GenePerc[i] = 100 - round(((temp$PauseStart[1] - temp$GeneStart[1]) / (temp$GeneEnd[1] - temp$GeneStart[1])) * 100, 2)
       data$DistFromTSS[i] = temp$GeneEnd[1] - temp$PauseStart[1]
       data$DistFrompA[i] = temp$PauseStart[1] - temp$GeneEnd[1]
    }


    # Melting temperature
    data$TmNN[i] = Tm_NN(data$Seq[i])

} # End of for loop

# Get DNA shape features (median value for 21 bp and value at nucleotide of pause)
# Read data
if (TYPE == "Real") {
   pred <- getShape(paste("../5_PauseLoci/", MUT, ".IDRrep.pauseWindows.fasta", sep = ""))
} else {
   pred <- getShape(paste("../5_PauseLoci/", MUT, ".IDRrep.pauseWindowsSHUFFLE.fasta", sep = ""))
}

data$ShapeEP = pred$EP[,11]
data$ShapeHelT = pred$HelT[,11]
data$ShapeMGW = pred$MGW[,11]
data$ShapeProT = pred$ProT[,11]
data$ShapeRoll = pred$Roll[,11]

# Add category label (real / shuffled)
data$Real = TYPE

# Report output
write.table(data, paste(MUT, ".IDRrep.PAUSESwithFEATURES_", TYPE, ".bed", sep=""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
