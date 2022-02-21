#!/n/app/R/4.0.1/bin/Rscript

### USE:
### sbatch -p short -t 0-00:30 -e logs/makeGOtables.err -o logs/makeGOtables.log --wrap="../Scripts/MakeGOtables.R"

library(data.table)

muts=c("bre1","bur2","bye1","cac1","cac2","cac3","cbc1","ccr4","cdc39","cdc73","chd1","ctr9","dhh1","dst1","eaf1","elf1","gcn5","hda1","hpc2","htz1","ino80","isw1","isw2","leo1","nap1","npl3","paf1","rad6","rco1","rpb4","rph1","rsc30","rtf1","rtt103","set2","set3","spt4","swr1","ubp8","vps15","vps34")


args<-commandArgs(TRUE)

method <- "Parent-Child-Intersection"

for (mut in muts) {
UPtable <- data.table(read.table(paste0('table-',mut,'.DEgeneList.UP-',method,'-None.txt'), header=TRUE))
DOWNtable <- data.table(read.table(paste0('table-',mut,'.DEgeneList.DOWN-',method,'-None.txt'),header=TRUE))
ALLtable <- data.table(read.table(paste0('table-',mut,'.DEgeneList.ALL-',method,'-None.txt'),header=TRUE))

# Calculate expected number of hits
UPtable[, Exp := round((Pop.term/Pop.total)*Study.total, 2)]
DOWNtable[, Exp := round((Pop.term/Pop.total)*Study.total, 2)]
ALLtable[, Exp := round((Pop.term/Pop.total)*Study.total, 2)]

# Calculate fold enrichment
UPtable[, FE := round(Study.term/Exp, 1)]
DOWNtable[, FE := round(Study.term/Exp, 1)]
ALLtable[, FE := round(Study.term/Exp, 1)]

# Filter
UPtableFilt <- UPtable[FE > 2 & p.adjusted < 0.05]
DOWNtableFilt <- DOWNtable[FE > 2 & p.adjusted < 0.05]
ALLtableFilt <- ALLtable[FE > 2 & p.adjusted < 0.05]

# Make new data tables
UP <- data.table(term = UPtableFilt$name, mut = paste0(UPtableFilt$FE,' (',round(UPtableFilt$p.adjusted, 6), ')'))
colnames(UP) <- c('term', mut)

DOWN <- data.table(term = DOWNtableFilt$name, mut = paste0(DOWNtableFilt$FE,' (',round(DOWNtableFilt$p.adjusted, 6), ')'))
colnames(DOWN) <- c('term', mut)

ALL <- data.table(term = ALLtableFilt$name, mut = paste0(ALLtableFilt$FE,' (',round(ALLtableFilt$p.adjusted, 6), ')'))
colnames(ALL) <- c('term', mut)

# Assign to table
assign(paste0(mut,'UP'), UP)
assign(paste0(mut,'DOWN'), DOWN)
assign(paste0(mut,'ALL'), ALL)

}

# Make list of DTs
UP_DTlist=list()
for (mut in muts) {
UP_DTlist <- append(UP_DTlist, eval(paste0(mut, 'UP')))
}
UP_DTs <- lapply(UP_DTlist, function(x) get(x))

DOWN_DTlist=list()
for (mut in muts) {
DOWN_DTlist <- append(DOWN_DTlist, eval(paste0(mut, 'DOWN')))
}
DOWN_DTs <- lapply(DOWN_DTlist, function(x) get(x))

ALL_DTlist=list()
for (mut in muts) {
ALL_DTlist <- append(ALL_DTlist, eval(paste0(mut, 'ALL')))
}
ALL_DTs <- lapply(ALL_DTlist, function(x) get(x))

# Merge by GO term
UP_all <- Reduce(
  function(x, y) merge(x, y, by='term', all=TRUE), 
  UP_DTs
)

DOWN_all <- Reduce(
  function(x, y) merge(x, y, by='term', all=TRUE), 
  DOWN_DTs
)

ALL_all <- Reduce(
  function(x, y) merge(x, y, by='term', all=TRUE), 
  ALL_DTs
)

# Write tables
write.table(UP_all, file='UP_GO_table.txt', row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
write.table(DOWN_all, file='DOWN_GO_table.txt', row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')
write.table(ALL_all, file='ALL_GO_table.txt', row.names=FALSE, col.names=TRUE, quote=FALSE, sep='\t')

