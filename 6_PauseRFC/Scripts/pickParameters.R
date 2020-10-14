#!/n/app/R/3.6.1/bin/Rscript

# Find optimal parameters for a random forest classifier and then predict pauses on test set

# Written by: K. Lachance
# Inspired by code: https://www.r-bloggers.com/how-to-implement-random-forests-in-r/
# Date: April 22, 2020

# Use: ./pickParameters.R mut nCat

# Load libraries
cat("Loading libraries...\n")
library(randomForest)

# Set seed for reproducibility
set.seed(123)

# Get input parameters
args <- commandArgs(trailingOnly = TRUE)
MUT <- args[1]
NCAT <- as.numeric(args[2]) # Number of categories in which to break continuous features

# Load data
cat("Formatting data...\n")
data_in = read.table(paste(MUT, ".IDRrep.PAUSESwithFEATURES.bed", sep=""), header = TRUE, stringsAsFactors = TRUE)

# Turn all entries into factors
data_factors = c('Real', 'PauseNucA', 'PauseNucB', 'PauseNucC', 'PauseNucD', 'PauseNucE', 'PauseNucAB', 'PauseNucBC', 'PauseNucCD', 'PauseNucDE')
data = data_in[,data_factors]

N = length(data_factors)
data_noFactors = data_in[,!names(data_in) %in% c(data_factors, 'Chr', 'Start', 'End', 'Strand', 'Seq')]

for (i in 1:ncol(data_noFactors)) {
    COL = N + i
    DATA_TYPE = colnames(data_noFactors)[i]
    tmp = data_noFactors[,i]
    tmp[is.na(tmp)] = 0
    data[, COL] = cut(tmp, NCAT) 
    colnames(data)[COL] <- paste(DATA_TYPE)
}

# Split into training (75%) and testing (25%) data
cat("Splitting training and testing data...\n")

### ADDED FOR SELECTING SPECIFIC FEATURES
#FOI = c("PauseNucA", "PauseNucE", "PauseNucAB", "PauseNucBC", "PauseNucCD", "PauseNucDE", "pCG", "h3k79me3", "h4k20me", "h3k36me3", "h3k79me", "h4r3me2s", "GenePerc", "TmNN", "ShapeHelT", "ShapeMGW", "ShapeProT", "ShapeRoll", "Ser5P") 
#DATA = data[,c("Real", FOI)]
#data = DATA

train <- sample(nrow(data), 0.75*nrow(data), replace = FALSE)
TrainSet <- data[train,]
TestSet <- data[-train,]

# Create a data frame to hold results
out <- data.frame(mtry = rep(c(1, 2, 3, 4, 5, 10, 15, 20, 25, 30), each = 4), ntree = rep(c(1000, 1500, 2000, 2500), 10), ErrMean = 0, ErrSD = 0, stringsAsFactors = FALSE)

# Create a Random Forest model with default parameters
cat("Trying different parameters...\n")
for (i in 1:nrow(out)) {

    # Get parameters 
    MTRY = out$mtry[i]
    NTREE = out$ntree[i]

    cat(paste("\tCreating random forest classifier with mtry = ", MTRY, " and ntree = ", NTREE, "\n", sep=""))

    model <- randomForest(Real ~ ., data = TrainSet, ntree = NTREE, mtry = MTRY, importance = TRUE)
    out$ErrMean[i] = mean(model$err.rate)
    out$ErrSD[i] = sd(model$err.rate)

} # End of for loop

# Write output to file
write.table(out, paste(MUT, ".randomForestParameters", NCAT, ".txt", sep=""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
