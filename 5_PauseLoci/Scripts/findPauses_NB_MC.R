#!/n/app/R/4.0.1/bin/Rscript

# Find pauses for given gene with specified NET-seq coverage
# Method defined in Churchman & Weissman (2011)

# USES NB DISTRIBUTION RATHER THAN GAUSSIAN

# Written by: K. Lachance
# Date: November 1, 2018
# Updated by: M. Couvillion 5/2021
# - Remove gene "buffer" at beginning and end that excludes those regions from pause calling
#   instead make window up/down size change if close to ends, to keep total window size constant


# Use: ./findPauses_NB_MC.R tmpCov.bed mut p w strand

# Import library
#install.packages("fitdistrplus", repos = "http://cran.us.r-project.org")
#install.packages("lsei", repos = "http://cran.us.r-project.org")
suppressMessages(library(fitdistrplus))
# suppressMessages(library(pscl))


# Read in parameters (data)
args <- commandArgs(trailingOnly = TRUE)
dat_file <- args[1]
mut <- args[2]
genecov <- args[3]
p <- as.numeric(args[4]) # SD over mean
w <- as.numeric(args[5]) # Half window size (nt on each side of potential pause)
strand <- args[6]

# Read file
dat <- read.table(dat_file)
colnames(dat) <- c("Chr", "Start", "End", "Cov")

# Function to calculate threshold for pauses given regional count data
calcThreshold <- function(x, p) {

  # Fit count data to negative binomial
  
#   if (min(x) == 0) {
#   
#   		fd <- zeroinfl(x~1|1,dist="negbin")
#      	size = 1
#    		mu=exp(coef(fd)[[1]])
#    		
#   } else {
  
  		fd = fitdist(x, "nbinom")
 		size = fd$estimate[[1]]
		mu = fd$estimate[[2]]
		
#   }
#   summary(zifd)
#   fit.coef <- coef(zifd)

    
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
# gene_length_buffer = nrow(dat)
# dat_noBuffer = dat[101:(gene_length_buffer-99),]

# potentialPauses = dat_noBuffer[dat_noBuffer$Cov > covT,]
potentialPauses = dat[dat$Cov > covT,]

# find start and end coords of gene
txptStart = dat$End[1]
txptEnd = dat$End[nrow(dat)]

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

	pp = potentialPauses$End[i]
	
    	# Extract in a +- w bp window surrounding position
    	if (pp > txptStart + w & pp < txptEnd - w) {
    	pauseRegionUp = dat[(dat$End > pp - w) & (dat$End < pp),]
    	pauseRegionDown = dat[(dat$End > pp) & (dat$End < pp + w),]
		} else if ( pp <= txptStart + w ) { 
		pauseRegionUp = dat[(dat$End < pp),]
		remainder = (2*w) - (pp - txptStart)
    	pauseRegionDown = dat[(dat$End > pp) & (dat$End < pp + remainder),]
    	} else if ( pp >= txptEnd - w ) {
		remainder = (2*w) - (txptEnd - pp)
		pauseRegionUp = dat[(dat$End < pp) & (dat$End > pp - remainder),]
    	pauseRegionDown = dat[(dat$End > pp),]
    	}
		
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

   # Write pauses to file
   write.table(pauses, file = paste(mut,"_cov",genecov, ".pausesNB_T", p, ".w", w, ".", strand, ".bed", sep=""), append = TRUE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

}