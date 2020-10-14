#!/n/app/R/3.6.1/bin/Rscript

# Find pauses for given gene with specified NET-seq coverage
# Method defined in Churchman & Weissman (2011)

# USES NB DISTRIBUTION RATHER THAN GAUSSIAN

# Written by: K. Lachance
# Date: November 1, 2018

# Use: ./findPauses_NB.R tmpCov.bed mut p strand

# Import library
#install.packages("fitdistrplus", repos = "http://cran.us.r-project.org")
#install.packages("lsei", repos = "http://cran.us.r-project.org")
suppressMessages(library(fitdistrplus))

# Read in parameters (data)
args <- commandArgs(trailingOnly = TRUE)
dat_file <- args[1]
mut <- args[2]
p <- as.numeric(args[3])
strand <- args[4]

# Read file
dat <- read.table(dat_file)
colnames(dat) <- c("Chr", "Start", "End", "Cov")

# Function to calculate threshold for pauses given regional count data
calcThreshold <- function(x, p) {

  # Fit count data to negative binomial
  fd = fitdist(x, "nbinom")
  
  # Get size and mu estimated parameters
  size = fd$estimate[[1]]
  mu = fd$estimate[[2]]
  
  # Calculate variance from size and mu
  var = mu + (mu^2) / size
  sd = sqrt(var)
  
  # Calculate threshold as mean + 3 * sd
  t = mu + p * sd
  return(t)
}

# Set minimum coverage threshold for a postion to potentially be considered a pause
covT = 2 # Changed to 2 (4 in original NET-seq paper)

# Create data frame containing all potential pauses (coverage > covT) within coding region (removing +- 100 bp buffer)
gene_length_buffer = nrow(dat)
dat_noBuffer = dat[101:(gene_length_buffer-99),]

potentialPauses = dat_noBuffer[dat_noBuffer$Cov > covT,]

# Create empty dataframe to hold pause sites
pauses = data.frame(Chr=character(), Start=integer(), End=integer(), Cov=integer(), stringsAsFactors = FALSE)

# While still finding pauses (looking via recursion)
lookingForPauses = TRUE

# If there are no potential pauses due to coverage theshold, do not look for pauses
if (nrow(potentialPauses) == 0) {
   lookingForPauses = FALSE
}

while (lookingForPauses) {

    # Set lookingForPauses to FALSE unless pause is found this round
    lookingForPauses = FALSE

    # For each potential pause
    for (i in 1:nrow(potentialPauses)) {

	pp = potentialPauses$Start[i]
	
    	# Extract in a +- 100 bp window surrounding position
    	pauseRegionUp = dat[(dat$Start > pp - 101) & (dat$Start < pp),]
    	pauseRegionDown = dat[(dat$Start > pp) & (dat$Start < pp + 101),]

    	# Remove any positions that are currently considered pauses (based on start position)
    	pauseRegionUp <- pauseRegionUp[ ! pauseRegionUp$Start %in% pauses$Start, ]
    	pauseRegionDown <- pauseRegionDown[ ! pauseRegionDown$Start %in% pauses$Start, ]

    	# Calculate threshold based on fit negative binomial distribution
    	pauseRegionCov = rbind(pauseRegionUp, pauseRegionDown)$Cov
    	t = calcThreshold(pauseRegionCov, p)

    	# If pp > t, add to list of pauses
    	# Otherwise, move to next possible pause
    	if (potentialPauses$Cov[i] > t) {
       	   pauses = rbind(pauses, potentialPauses[i,])    
	   lookingForPauses = TRUE
    	}
    }

    # Remove pauses from potentialPauses
    potentialPauses <- potentialPauses[ ! potentialPauses$Start %in% pauses$Start, ]
    
    # Exit while loop if no more potential pauses
    if (nrow(potentialPauses) == 0) {
       lookingForPauses = FALSE
    }

}

# If pauses were detected, sort and report them
if (nrow(pauses) > 0) {

   # Sort pauses by start position
   pauses <- pauses[order(pauses$Start),] 

   # Write pauses to standard output
   write.table(pauses, file = paste(mut, ".pausesNB_T", p, ".", strand, ".bed", sep=""), append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

}