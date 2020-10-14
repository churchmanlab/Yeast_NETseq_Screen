#!/n/app/R/3.6.1/bin/Rscript

# Shuffle pauses, maintaining the same number per gene, none overlapping each other or real pauses

# Written by: K. Lachance
# Date: May 6, 2020

# Use: ./shufflePauses.R mut strand

# Read in parameters (data)
args <- commandArgs(trailingOnly = TRUE)
MUT <- args[1]
STRAND <- args[2]

# Convert strand to symbol
if (STRAND == "pos") {
   strandSymbol = "+"
} else {
   strandSymbol = "-"
}

# Command to extract data in each window for each pause
bashCommand='/n/app/bedops/2.4.30/bedextract $1 $2 '

# Get data
pauseCounts = read.table(paste(MUT, ".pausesByGene.", STRAND, ".bed", sep=""), stringsAsFactors = FALSE)
colnames(pauseCounts) = c("Chr", "Start", "End", "N")

# Remove rows with 0 pauses in gene
pauseCounts = pauseCounts[pauseCounts$N > 0,]

# Create data frame to hold output
out = data.frame(Chr = character(), Start = integer(), End = integer(), Strand = character(), stringsAsFactors = FALSE)

# Loop through each gene
for (i in 1:nrow(pauseCounts)) {

    # Counter
    if (i %% 100 == 0) {cat(paste("Done with ", i, " / ", nrow(pauseCounts), " genes!\n"))}

    # Get number of pauses 
    N = pauseCounts$N[i]

    # Get real pauses in gene
    # Write temporary file for overlaps
    tmp = data.frame(Chr = pauseCounts$Chr[i], Start = pauseCounts$Start[i], End = pauseCounts$End[i], stringsAsFactors = FALSE)
    write.table(tmp, paste(MUT, ".tmpCoord.txt", sep=""), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

    realP = as.vector(read.table(pipe(paste(bashCommand, ' ', MUT, '.IDRrep.PAUSES.", STRAND, ".bed ',  MUT, '.tmpCoord.txt', sep='')), sep="\t", header=FALSE, col.names=c("Chr", "Start", "End", "Strand"), colClasses=c('character', 'numeric', 'numeric', 'character'))$Start)

    # Pick shuffled pauses
    nPausesPicked = 0
    nPausesLeft = N
    shuffleP = c()
    while (nPausesLeft > 0) {
        tempShuffleP = sample(pauseCounts$Start[i]:pauseCounts$End[i], nPausesLeft)
	shuffleP = c(shuffleP, tempShuffleP[!(tempShuffleP %in% realP)])
	nPausesPicked = length(shuffleP)
	nPausesLeft = N - nPausesPicked
	
    } # End of while loop

    tmpOut = data.frame(Chr = rep(pauseCounts$Chr[i], N), Start = sort(shuffleP), End = sort(shuffleP) + 1, Strand = rep(strandSymbol, N), stringsAsFactors = FALSE)
    out = rbind(out, tmpOut)

} # End of for loop
write.table(out, paste(MUT, ".IDRrep.SHUFFLE.", STRAND, ".bed", sep=""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)