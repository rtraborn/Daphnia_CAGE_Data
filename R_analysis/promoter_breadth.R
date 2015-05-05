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

dmel_widths <- width(dmel.gr)

widths.dpul <- as.data.frame(dpul_widths)
widths.dmel <- as.data.frame(dmel_widths)

names(widths.dpul) <- c("widths")
names(widths.dmel) <- c("widths")

#plotting the figures

f <- ggplot(widths.dpul, aes(widths))
f + geom_density()
f <- ggplot(widths.dmel, aes(widths))
f + geom_density()

ggsave("width_compare.png",width=6,height=5)

q()


