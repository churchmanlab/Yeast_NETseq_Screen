#!/n/app/R/3.6.1/bin/Rscript

# Written by: K. Lachance
# Date: June 9, 2020

# Make a CDF plot for all PIs for wildtype

# Use: ./plotCDFpi.R

library(ggplot2)
library(svglite)

data1 = read.table("wt_TSS.PI.bed", stringsAsFactors = F)[,c(6,13,14)]
colnames(data1) = c("Gene", "PI", "Region")

data2 = read.table("wt_pA.PI.bed", stringsAsFactors = F)[,c(6,13,14)]
colnames(data2) = c("Gene", "PI", "Region")

data3 = read.table("wt_5pSS.PI.bed", stringsAsFactors = F)[,c(6,13,14)]
colnames(data3) = c("Gene", "PI", "Region")

data4 = read.table("wt_3pSS.PI.bed", stringsAsFactors = F)[,c(6,13,14)]
colnames(data4) = c("Gene", "PI", "Region")

data5 = read.table("wt_Antisense.PI.bed", stringsAsFactors = F)[,c(6,13,14)]
colnames(data5) = c("Gene", "PI", "Region")

data = rbind(data1, data2, data3, data4, data5)


p = ggplot(data, aes(PI, group = Region, color = Region)) + 
  stat_ecdf(geom = "step") +
  scale_color_manual(values = c("TSS" = "#006838", "pA" = "#BE1E2D", "5pSS" = "#2E3192", "3pSS" = "#00AEEF", "Anti" = "#517BA8")) + 
  theme_bw() + 
  xlim(0, 5) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none")

ggsave(filename = "PI_CDF.svg", plot = p, height = 5, width = 5)

