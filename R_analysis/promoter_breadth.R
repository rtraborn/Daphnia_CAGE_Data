###calculates promoter breadth in Daphnia and Drosophila using CAGE annotations

require("GenomicRanges")
require("ggplot2")

dpul.df <- read.table(file="/home/rtraborn/research/Daphnia/Daphnia_CAGE_Data/promoter_arch/D_pulex/All.samples.tagClusters.qLow0.1_qUp0.9.bed",header=F,stringsAsFactors=F)

names(dpul.df) <- c("chr","start","end","score","score2","strand","start2","end2","score3","score4","score5","score6")
                                                                                                                                      
dmel.df <- read.table(file="/home/rtraborn/research/Daphnia/Daphnia_CAGE_Data/promoter_arch/dmel_hoskins/All.samples.tagClusters.qLow0.1_qUp0.9.bed",header=F,stringsAsFactors=F)

names(dmel.df) <- c("chr","start","end","score","score2","strand","start2","end2","score3","score4","score5","score6")

dpul.gr <- makeGRangesFromDataFrame(dpul.df, seqnames.field="chr", start.field="start", end.field="end", strand.field="strand")

dmel.gr <- makeGRangesFromDataFrame(dmel.df, seqnames.field="chr", start.field="start", end.field="end", strand.field="strand")

dpul_widths <- width(dpul.gr)
n_widths_dp <- length(dpul_widths)

dmel_widths <- width(dmel.gr)
n_widths_dm <- length(dmel_widths)

n_diff <- n_widths_dp-n_widths_dm 

print(length(dpul_widths))

dmel_widths2 <-c(dmel_widths,rep(NA,n_diff))

print(length(dmel_widths2))

widths.df  <- as.data.frame(cbind(dpul_widths,dmel_widths2))
names(widths.df) <- c("D_pulex","D_melanogaster")

#plotting the figures

f <- ggplot(widths.df, aes(dpul_widths, dmel_widths2))
f + geom_boxplot()


ggsave("width_compare.png",width=6,height=5)

q()


