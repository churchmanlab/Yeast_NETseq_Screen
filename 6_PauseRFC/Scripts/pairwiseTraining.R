#!/n/app/R/3.6.1/bin/Rscript

# Train and predict pauses pairwise across deletion strains

# Written by: K. Lachance
# Date: May 7, 2020

# Use: ./pairwiseTraining.R deletionStrains.txt MTRY NTREE

# Load libraries
cat("Loading libraries...\n")
library(randomForest)
library(ROCR)
library(dplyr)

# Set seed for reproducibility
set.seed(123)

# Get input parameters
args <- commandArgs(trailingOnly = TRUE)
m <- args[1]
muts2 <- as.vector(read.table(args[2], stringsAsFactors = FALSE)$V1)
MTRY <- as.numeric(args[3])
NTREE <- as.numeric(args[4])

# Double loop for each combination of deletion strain
#for (m in muts1) {
    
    for (n in muts2) {

    	if (m == n) {next}

        # Report pair
    	cat(paste("Working on ", m,	" - ", n, " pair!\n", sep=""))

    	# Load data
    	data1_in = read.table(paste(m, ".IDRrep.PAUSESwithFEATURES.bed", sep=""), header = TRUE, stringsAsFactors = TRUE)
    	data2_in = read.table(paste(n, ".IDRrep.PAUSESwithFEATURES.bed", sep=""), header = TRUE, stringsAsFactors = TRUE)
    
	# Remove pauses in testing data set present in training data
	data2_temp = anti_join(data2_in, data1_in, by=c('Chr', 'Start', 'End', 'Start'))
	data2_in = data2_temp

	# Turn all entries into factors
    	data_factors = c('Real', 'PauseNucA', 'PauseNucB', 'PauseNucC', 'PauseNucD', 'PauseNucE', 'PauseNucAB', 'PauseNucBC', 'PauseNucCD', 'PauseNucDE')
    	data1 = data1_in[,data_factors]
    	data2 = data2_in[,data_factors]

    	N = length(data_factors)
    
	data1_noFactors = data1_in[,!names(data1_in) %in% c(data_factors, 'Chr', 'Start', 'End', 'Strand', 'Seq')]
	data2_noFactors = data2_in[,!names(data2_in) %in% c(data_factors, 'Chr', 'Start', 'End', 'Strand', 'Seq')]
	
    	for (i in 1:ncol(data1_noFactors)) {
    	    COL = N + i
    	    DATA_TYPE = colnames(data1_noFactors)[i]
    	    tmp = c(data1_noFactors[,i], data2_noFactors[,i])
    	    tmp[is.na(tmp)] = 0
    	    data1[, COL] = cut(tmp, 5, labels = c("a", "b", "c", "d", "e"))[1:nrow(data1)]
	    data2[, COL] = cut(tmp, 5, labels = c("a", "b", "c", "d", "e"))[(nrow(data1)+1):length(tmp)]
    	    colnames(data1)[COL] <- paste(DATA_TYPE)
	    colnames(data2)[COL] <- paste(DATA_TYPE)
    	}

	### ADDED FOR SELECTING SPECIFIC FEATURES
	#FOI = c("PauseNucA", "PauseNucE", "PauseNucAB", "PauseNucBC", "PauseNucCD", "PauseNucDE", "pCG", "h3k79me3", "h4k20me", "h3k36me3", "h3k79me", "h4r3me2s", "GenePerc", "DistFromTSS", "TmNN", "ShapeHelT", "ShapeMGW", "ShapeProT", "ShapeRoll") #Old
#	FOI = c("PauseNucA", "PauseNucE", "PauseNucAB", "PauseNucBC", "PauseNucCD", "PauseNucDE", "pCG", "h3k79me3", "h4k20me", "h3k36me3", "h3k79me", "h4r3me2s", "GenePerc", "TmNN", "ShapeHelT", "ShapeMGW", "ShapeProT", "ShapeRoll", "Ser5P") # New
#	DATA1 = data1[,c("Real", FOI)]
#	DATA2 = data2[,c("Real", FOI)]
	data1 = DATA1
	data2 = DATA2

	# Train random forest classifier
	model <- randomForest(Real ~ ., data = data1, ntree = NTREE, mtry = MTRY, importance = TRUE)
	prediction_for_roc_curve <- predict(model,data2[,-1],type="prob")

	# Define which observations belong to class[i]
	true_values <- ifelse(data2[,1]=='Real',1,0)

	# Assess the performance of classifier for Real pauses
	pred <- prediction(prediction_for_roc_curve[,1],true_values)
	auc <- performance(pred, "auc")
	tmp2 = data.frame(Train = paste(m), Test = paste(n), AUC = unlist(slot(auc, "y.values")), stringsAsFactors = FALSE)
	write.table(tmp2, "AUCvalues_pairwise.txt", append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

    } # End of inside for loop
#} # End of outside for loop