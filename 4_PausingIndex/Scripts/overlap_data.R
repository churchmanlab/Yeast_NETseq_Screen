#!/n/app/R/3.6.1/bin/Rscript

# Overlap data (and fill in 0s) to make a matrix of data; for each row, include the NET-seq value at every nucleotide position

# Use: ./overlap_data.R data.bed regions.of.interest.bed strand
# Example: ./overlap_data.R wildtype.bedgraph PolIIGenes.Trimmed250.Plus.bed pos wildtype

# Load libraries
#install.packages("zoo", repos = "http://cran.us.r-project.org")
library(zoo)

# Allow arguments
args <- commandArgs(trailingOnly = TRUE)
dat_file <- args[1]
window_file <- args[2]
ST <- args[3]
mut <- args[4]

# Read tables
window <- read.table(window_file, colClasses=c("character", "integer", "integer", "character", "character", "character"))
colnames(window) = c("Chr", "Start", "End", "Gene", "Value", "Strand")
n = nrow(window)

# Calculate window size
w = window$End[1] - window$Start[1]

# Command to extract NETseq data in each window
bashCommand='/n/app/bedops/2.4.30/bedextract $1 $2 '

# Create matrix to hold expanded window output
dat_mt = data.frame(matrix(nrow=n, ncol=(w+4), 0), stringsAsFactors=FALSE)

# For each window
for (i in 1:n){

    # Get chromosome, window start, and window end
    chr = as.character(window$Chr[i])
    start = as.integer(window$Start[i])
    end = as.integer(window$End[i])
    gene = as.character(window$Gene[i])

    # Write temporary file
    tmp = data.frame(Chr=chr, Start=start, End=end)
    write.table(tmp, paste(mut, '_tmpCoord_', w, '.txt', sep=''), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

    # Call bash command
    temp = read.table(pipe(paste(bashCommand, dat_file, ' ', mut, '_tmpCoord_', w, '.txt', sep='')), sep="\t", header=FALSE, col.names=c("Chr", "Start", "End", "Value"), colClasses=c('character', rep('numeric', 3)))

    n_temp = dim(temp)[1]

    # Fill in unlisted positions in temp
    t = data.frame(Chr=character(), Start=numeric(), End=numeric(), Value=numeric(), stringsAsFactors=FALSE)

    # If there are data points in the bed files within the window of current interest
    if ( n_temp > 1 ) {
        # Fill in missing data points so that every base in the window is accounted for
        for (k in 1:(n_temp-1)) {
            if (temp[k,3] != temp[(k+1),2]) {
               l = c(temp[k, 1], temp[k, 3], temp[k+1, 2], 0)
               t[(nrow(t)+1),] = temp[k,]
               t[(nrow(t)+1),] = l
            }
            else {
               t[(nrow(t)+1),] = temp[k,]
            }
        }
        t[(nrow(t)+1),] = temp[k+1,]
        # Set t equal to temp
        temp <- t
    }

    # Reset start and ends of windows to be the actual start and end of window
    if (nrow(temp) > 0 ) {
        temp[1,'Start'] = start
        temp[nrow(temp),'End'] = end

        # Create vector of values for every position in window
        nbStarts = rep(as.numeric(temp[,'Value']), (as.numeric(temp[,'End']) - as.numeric(temp[,'Start'])))

        dat_mt[i, 1:4] <- c(chr, start, end, gene)

        # Reverse if on negative strand
        if (ST=='pos') {
            dat_mt[i,5:(w+4)] = nbStarts
        } else {
	     dat_mt[i,5:(w+4)] = rev(nbStarts)
        }
    }
}

write.table(dat_mt, file=paste(mut, "_OverlapMatrix.", w/2, ".", ST, ".txt", sep=""), sep="\t", col.names=FALSE, quote=FALSE, row.names=FALSE)
