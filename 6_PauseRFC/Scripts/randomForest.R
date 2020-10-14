#!/n/app/R/3.6.1/bin/Rscript

# Implement a random forest classifier and then predict pauses on test set

# Written by: K. Lachance
# Inspired by code: https://www.r-bloggers.com/how-to-implement-random-forests-in-r/
# Date: April 22, 2020

# Use: ./randomForest.R sample mtry ntree

# Load libraries
cat("Loading libraries...\n")
library(randomForest)
library(ROCR)
library(ggplot2)
library(svglite)

# Set seed for reproducibility
set.seed(123)

# Get input parameters
args <- commandArgs(trailingOnly = TRUE)
MUT <- args[1]
MTRY <- as.numeric(args[2])
NTREE <- as.numeric(args[3])

cat(paste("\t***", MUT, "***\n", sep=""))

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
    data[, COL] = cut(tmp, 5)
    colnames(data)[COL] <- paste(DATA_TYPE)
}

# Read in feature categories
featureCat = read.table("featureCategories.txt", stringsAsFactors = FALSE) 
colnames(featureCat) = c("Feat", "Cat")

# Create data frame to hold values used to plot ROC
outROC = data.frame(X = as.numeric(), Y = as.numeric(), LAB = as.character(), stringsAsFactors = FALSE)

# Loop through all categories
for (CAT in unique(featureCat$Cat)) {

   cat(paste("Collecting all ", CAT, " features...\n", sep=""))
   FEAT = featureCat[featureCat$Cat == CAT, 'Feat']

   DATA = data[,c("Real", FEAT)]

   # Split into training (75%) and testing (25%) data
   cat("\tSplitting training and testing data...\n")
   train <- sample(nrow(DATA), 0.75*nrow(DATA), replace = FALSE)
   TrainSet <- DATA[train,]
   TestSet <- DATA[-train,]

   cat("\tTraining random forest classifier...\n")
   model <- randomForest(Real ~ ., data = TrainSet, ntree = NTREE, mtry = MTRY, importance = TRUE)
   
   IMP = importance(model)
   write.table(IMP, paste(MUT, ".", CAT, "_importance.txt", sep=""), quote = FALSE, sep="\t")

   cat("\tPlotting ROC curves and calculating AUC...\n")

   prediction_for_roc_curve <- predict(model,TestSet[,-1],type="prob")

   # Define which observations belong to class[i]
   true_values <- ifelse(TestSet[,1]=='Real',1,0)

   # Assess the performance of classifier for Real pauses
   pred <- prediction(prediction_for_roc_curve[,1],true_values)
   perf <- performance(pred, "tpr", "fpr")
   tmp = data.frame(X = unlist(slot(perf, "x.values")), Y = unlist(slot(perf, "y.values")), LAB = paste(CAT), stringsAsFactors = FALSE)
   outROC = rbind(outROC, tmp)
  
   auc <- performance(pred, "auc")
   tmp2 = data.frame(FeatureCategory = paste(CAT), Mut = paste(MUT), AUC = unlist(slot(auc, "y.values")), stringsAsFactors = FALSE)
   write.table(tmp2, "AUCvalues.txt", append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

}

# Try combination of nucleotide sequence and chip-seq
cat(paste("Collecting all Seq + ChIP features...\n", sep=""))
FEAT = featureCat[featureCat$Cat == "Seq" | featureCat$Cat == "ChIP", 'Feat']
DATA = data[,c("Real", FEAT)]

# Split into training (75%) and testing (25%) data
cat("\tSplitting training and testing data...\n")
train <- sample(nrow(DATA), 0.75*nrow(DATA), replace = FALSE)
TrainSet <- DATA[train,]
TestSet <- DATA[-train,]

cat("\tTraining random forest classifier...\n")
model <- randomForest(Real ~ ., data = TrainSet, ntree = NTREE, mtry = MTRY, importance = TRUE)

IMP = importance(model)
write.table(IMP, paste(MUT, ".Seq+ChIP_importance.txt", sep=""), quote = FALSE, sep="\t")

cat("\tPlotting ROC curves and calculating AUC...\n")
prediction_for_roc_curve <- predict(model,TestSet[,-1],type="prob")

# Define which observations belong to class[i]
true_values <- ifelse(TestSet[,1]=='Real',1,0)

# Assess the performance of classifier for Real pauses
pred <- prediction(prediction_for_roc_curve[,1],true_values)
perf <- performance(pred, "tpr", "fpr")
tmp = data.frame(X = unlist(slot(perf, "x.values")), Y = unlist(slot(perf, "y.values")), LAB = "xxxCOMB", stringsAsFactors = FALSE)
outROC = rbind(outROC, tmp)

auc <- performance(pred, "auc")
tmp2 = data.frame(FeatureCategory = "Seq+ChIP", Mut = paste(MUT), AUC = unlist(slot(auc, "y.values")), stringsAsFactors = FALSE)
write.table(tmp2, "AUCvalues.txt", append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Try combination of nucleotide sequence and shape features
cat(paste("Collecting all Seq + Shape features...\n", sep=""))
FEAT = featureCat[featureCat$Cat == "Seq" | featureCat$Cat == "Shape", 'Feat']
DATA = data[,c("Real", FEAT)]

# Split into training (75%) and testing (25%) data
cat("\tSplitting training and testing data...\n")
train <- sample(nrow(DATA), 0.75*nrow(DATA), replace = FALSE)
TrainSet <- DATA[train,]
TestSet <- DATA[-train,]

cat("\tTraining random forest classifier...\n")
model <- randomForest(Real ~ ., data = TrainSet, ntree = NTREE, mtry = MTRY, importance = TRUE)

IMP = importance(model)
write.table(IMP, paste(MUT, ".Seq+Shape_importance.txt", sep=""), quote = FALSE, sep="\t")

cat("\tPlotting ROC curves and calculating AUC...\n")
prediction_for_roc_curve <- predict(model,TestSet[,-1],type="prob")

# Define which observations belong to class[i]
true_values <- ifelse(TestSet[,1]=='Real',1,0)

# Assess the performance of classifier for Real pauses
pred <- prediction(prediction_for_roc_curve[,1],true_values)
perf <- performance(pred, "tpr", "fpr")
tmp = data.frame(X = unlist(slot(perf, "x.values")), Y = unlist(slot(perf, "y.values")), LAB = "yyyCOMB", stringsAsFactors = FALSE)
outROC = rbind(outROC, tmp)

auc <- performance(pred, "auc")
tmp2 = data.frame(FeatureCategory = "Seq+Shape", Mut = paste(MUT), AUC = unlist(slot(auc, "y.values")), stringsAsFactors = FALSE)
write.table(tmp2, "AUCvalues.txt", append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Finally, combine all features
cat("Collecting ALL features...\n")

# Split into training (75%) and testing (25%) data
cat("\tSplitting training and testing data...\n")
train <- sample(nrow(data), 0.75*nrow(data), replace = FALSE)
TrainSet <- data[train,]
TestSet <- data[-train,]

cat("\tTraining random forest classifier...\n")
model <- randomForest(Real ~ ., data = TrainSet, ntree = NTREE, mtry = MTRY, importance = TRUE)

cat("\tPlotting ROC curves and calculating AUC...\n")

prediction_for_roc_curve <- predict(model,TestSet[,-1],type="prob")

# Define which observations belong to class[i]
true_values <- ifelse(TestSet[,1]=='Real',1,0)
# Assess the performance of classifier for Real pauses
pred <- prediction(prediction_for_roc_curve[,1],true_values)
perf <- performance(pred, "tpr", "fpr")
tmp = data.frame(X = unlist(slot(perf, "x.values")), Y = unlist(slot(perf, "y.values")), LAB = "zzzCOMB", stringsAsFactors = FALSE)
outROC = rbind(outROC, tmp)

auc <- performance(pred, "auc")
tmp2 = data.frame(FeatureCategory = "COMB", Mut = paste(MUT), AUC = unlist(slot(auc, "y.values")), stringsAsFactors = FALSE)
write.table(tmp2, "AUCvalues.txt", append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

# Plot ROC curve
cat("Plotting ROC curve...\n")

p <- ggplot(outROC, aes(x = X, y = Y, color = LAB)) + 
  geom_abline(intercept = 0, slope = 1, color = "grey50") + 
  geom_line() + 
  theme_bw() + 
  xlab("False positive rate") + ylab("True positive rate") + 
  scale_color_discrete(name = "Features Included", labels = c("Chromatin landscape", "PolII CTD Modifications", "Postion within Gene", "Nucleotide sequence", "DNA shape", "Chromatin + Sequence", "Shape + Sequence", "All features")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  ggtitle(paste(MUT, " ROC curve\nmtry = ", MTRY, ", ntrees = ", NTREE, sep=""))

ggsave(filename = paste(MUT, ".ROC.svg", sep=""), plot = p, width = 6, height = 4)


# Outputting confusion matrix and importance of each feature
cat("Calculating model accuracy and importance of each feature...\n")
pred = predict(model, newdata=TestSet[,-1])
cm = table(TestSet[,1], pred)
print(cm)

IMP = importance(model)
write.table(IMP, paste(MUT, ".importance.txt", sep=""), quote = FALSE, sep = "\t", row.names = TRUE)

cat("Done!\n\n\n")
