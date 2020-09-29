#!/n/app/R/3.6.1/bin/Rscript

# Generate boxplot for AS/S ratio for all mutants

# Written by: K. Lachance
# Date: October 30, 2018

# Use: ./plotAStoSboxplot.R ALL.NETseqRatio.bed coding/tandem

# Load libraries
library(ggplot2)
library(reshape2)
library(svglite)

# Read in parameters (data)
args <- commandArgs(trailingOnly = TRUE)
dat_file <- args[1]
label <- args[2]

# Read file
dat <- read.table(dat_file, stringsAsFactors = FALSE)
colnames(dat) <- c("Chr", "Start", "End", "Gene", "Score", "Strand", "Sense", "Antisense", "Ratio", "Mut")

# Remove rows in which there is no sense expression (Ratio = NA) or antisense expression (Ratio = 0)
dat <- dat[!(is.na(dat$Ratio)),]
dat <- dat[(dat$Ratio != 0),]

# Calculate the log2(Ratio)
dat$LogRatio = log2(dat$Ratio)

# Calculate median AS/S ratio for each strain
medians = aggregate(dat$LogRatio, list(dat$Mut), median, rm.na = TRUE)
colnames(medians) = c("Mut", "Med")

# Report medians
if (label == "coding") {
   filename = "CodingAntisense/AntisenseSenseRatioMedians.txt"
} else {
   filename = "TandemAntisense/AntisenseSenseRatioMedians.txt"
}
write.table(medians, paste(filename), quote = FALSE, col.names = FALSE, row.names = FALSE, sep="\t")


# Reorder factors by median (small > large)
strain_order = medians[order(medians$Med),'Mut']
dat$Mut = ordered(dat$Mut, levels = strain_order)

# Add coloring label 
wt_median = medians[(medians$Mut == 'wt'), 'Med']
wt_lower = quantile(dat[(dat$Mut == 'wt'), 'LogRatio'], 0.45)
wt_upper = quantile(dat[(dat$Mut == 'wt'), 'LogRatio'], 0.55)

dat$Lab = "none" # Default color label is 0 (only for wildtype)
dat[(dat$Mut == 'wt'), 'Lab'] = "wt"

for (mut in strain_order) { # For each mut
    
    if (medians[(medians$Mut == mut), 'Med'] < wt_lower) { # If mut median is less than wildtype
	dat[(dat$Mut == mut), 'Lab'] = "S" # Label with 1 (more antisense)
    } else if (medians[(medians$Mut ==	mut), 'Med'] > wt_upper) { # If mut median is greater than wildtyle
      	dat[(dat$Mut == mut), 'Lab'] = "AS" # Label with -1 (less antisense)
    }
 
}

# Calculate standard deviation of wildtype for plotting
wt_sd = sd(dat[(dat$Mut == 'wt'), 'LogRatio'])

# Plot
if (label == "coding") {
   filename = "boxplotAStoSratio_coding.svg"
   plotTitle = "Gene Body Antisense Transcription"
} else {
   filename = "boxplotAStoSratio_tandem.svg"
   plotTitle = "Divergent Antisense Transcription"
}

p1 <- ggplot(dat, aes(x = Mut, y = LogRatio, group = Mut)) + 
   geom_boxplot(outlier.shape = NA, coef = 0, notch = TRUE, aes(fill = Lab)) + # Added coef to not plot whiskers
   scale_fill_manual(values=c("#507BA7", "#FFFFFF", "#A7533F", "#6D6E71")) + # For coding & tandem
   geom_hline(aes(yintercept=wt_median), color="black") +
   geom_hline(aes(yintercept=wt_lower), color="black", linetype = "dotted") +
   geom_hline(aes(yintercept=wt_upper), color="black", linetype = "dotted") +
   coord_cartesian(ylim = c(-6, 0)) + # Adjust as necessary
   theme_bw() + 
   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
   theme(legend.position="none") + 
   ylab("log2(Antisense/Sense Ratio)") +
   ggtitle(paste(plotTitle)) 

ggsave(file=paste(filename), plot=p1, width=10, height=6)

