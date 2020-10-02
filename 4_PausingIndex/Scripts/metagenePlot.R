#!/n/app/R/3.6.1/bin/Rscript

# Written by: K. Lachance
# Date: August 17, 2017

# Plot the expanded bed files around loci of interest

# Use: ./metagenePlot.R mut

# Import libraries
library(ggplot2)
library(svglite)

# Get arguments (matrices)
args <- commandArgs(trailingOnly = TRUE)
mut <- args[1]

# Combine data by deletion string
mut_data <- read.table(paste("./readyToPlot/", mut, "_combinedMeltedMatrix.bed", sep=""), stringsAsFactors = FALSE)
wt_data <- read.table("./readyToPlot/wt_combinedMeltedMatrix.bed", stringsAsFactors = FALSE)

data = rbind(mut_data, wt_data)
colnames(data) = c("Gene", "Pos", "Value", "Mut", "Region")

# Move wt to beginning (alphabetically)
data[data$Mut == "wt", 'Mut'] = "aaaWT"

# Plot
p_tss <- ggplot(data[data$Region == "TSS",], aes(Pos, Value*1000, color = Mut, fill = Mut)) +
   geom_smooth(span = 0.01, na.rm = TRUE) +
   scale_color_manual(values=c("grey75", "#006838")) + scale_fill_manual(values=c("grey75", "#006838")) + 
   theme_bw() +
   geom_vline(xintercept = 100, linetype="longdash") +
   scale_x_continuous(breaks = c(100, 350, 600), labels = c(0, 250, 500)) +
   coord_cartesian(ylim =c(0.5, 3)) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
   theme(legend.position = "none")

ggsave(file=paste(mut, "vWT_metagene_TSS.svg", sep=""), plot=p_tss, width=5, height=3)

p_pA <- ggplot(data[data$Region == "pA",], aes(Pos, Value*1000, color = Mut, fill = Mut)) +
   geom_smooth(span = 0.01, na.rm = TRUE) +
   scale_color_manual(values=c("grey75", "#be1e2d")) + scale_fill_manual(values=c("grey75", "#be1e2d")) +
   theme_bw() +
   geom_vline(xintercept = 500, linetype="longdash") +
   scale_x_continuous(breaks = c(0, 250, 500)) +
   coord_cartesian(ylim =c(1, 2.5)) +
   scale_y_continuous(breaks = c(1, 2)) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
   theme(legend.position = "none")

ggsave(file=paste(mut, "vWT_metagene_pA.svg", sep=""), plot=p_pA, width=5, height=3)

p_anti <- ggplot(data[data$Region == "Anti",], aes(Pos, Value*1000, color = Mut, fill = Mut)) +
   geom_smooth(span = 0.01, na.rm = TRUE) +
   scale_color_manual(values=c("grey75", "#507BA6")) + scale_fill_manual(values=c("grey75", "#507BA6")) +
   theme_bw() +
   geom_vline(xintercept = 500, linetype="longdash") +
   scale_x_continuous(breaks = c(0, 250, 500)) +
   coord_cartesian(ylim =c(0.5, 2.5)) +
   scale_y_continuous(breaks = c(1,2)) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
   theme(legend.position = "none")

ggsave(file=paste(mut, "vWT_metagene_Anti.svg", sep=""), plot=p_anti, width=5, height=3)

# Split data by region
SS5p_data = data[data$Region == "5pSS",]
SS3p_data = data[data$Region == "3pSS",]

SS3p_data$Pos = SS3p_data$Pos + 55

SS5p_data$Cat = "NA"
SS3p_data$Cat = "NA"

SS5p_data[SS5p_data$Mut == "aaaWT", 'Cat'] = "aaa"
SS5p_data[SS5p_data$Mut != "aaaWT", 'Cat'] = "bbb"
SS3p_data[SS3p_data$Mut == "aaaWT", 'Cat'] = "ccc"
SS3p_data[SS3p_data$Mut != "aaaWT", 'Cat'] = "ddd"

SS_data = rbind(SS5p_data, SS3p_data)

# Plot
p_SS <- ggplot(SS_data, aes(Pos, Value*100, color = Cat, fill = Cat, group = Cat)) +
   geom_smooth(span = 0.01) +
   scale_color_manual(values=c("grey75", "#2e3192", "grey75", "#00aeef")) + scale_fill_manual(values=c("grey75", "#2e3192","grey75", "#00aeef")) +
   theme_bw() +
   geom_vline(xintercept = c(25, 80), linetype="longdash") +
   scale_x_continuous(breaks = c(25, 80)) +
   coord_cartesian(ylim =c(0, 6)) +
   scale_y_continuous(breaks = c(2,4,6)) +
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
   theme(legend.position = "none") 

ggsave(file=paste(mut, "vWT_metagene_SS.svg", sep=""), plot=p_SS, width=5, height=3)

