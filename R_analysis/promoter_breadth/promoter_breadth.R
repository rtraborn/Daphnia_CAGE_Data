###calculates promoter breadth in Daphnia and Drosophila using CAGE annotations

require("GenomicRanges")
require("ggplot2")
require("scales")

dpul.df <- read.table(file="/home/rtraborn/Daphnia/Daphnia_CAGE_Data/promoter_arch/D_pulex/ConsensusClusters_PA42.bed",header=F,stringsAsFactors=F,skip=1)

names(dpul.df) <- c("chr","start","end","placeholder","score","strand","start_core","end_core","score2")
                                                                                                                                      
dmel.df <- read.table(file="/home/rtraborn/Daphnia/Daphnia_CAGE_Data/promoter_arch/dmel_hoskins/All.samples.tagClusters.qLow0.1_qUp0.9.bed",header=F,stringsAsFactors=F)

names(dmel.df) <- c("chr","start","end","score","score2","strand","start2","end2","score3","score4","score5","score6")

dpul.gr <- makeGRangesFromDataFrame(dpul.df, seqnames.field="chr", start.field="start", end.field="end", strand.field="strand")

dmel.gr <- makeGRangesFromDataFrame(dmel.df, seqnames.field="chr", start.field="start", end.field="end", strand.field="strand")

dpul_widths <- width(dpul.gr)
n_widths_dp <- length(dpul_widths)

dmel_widths <- width(dmel.gr)
n_widths_dm <- length(dmel_widths)

width_vec <- c(dpul_widths, dmel_widths)
org_vec <- c(rep("D_pulex",n_widths_dp),rep("D_melanogaster",n_widths_dm))

widths.df <- data.frame(width_vec, org_vec)
names(widths.df) <- c("promoter_breadth","species")

#plotting the figures

ggplot(widths.df, aes(x=promoter_breadth,y = (..count..)/sum(..count..))) + geom_histogram(binwidth=5,colour="blue",fill="white") + facet_grid(species ~ .) + scale_y_continuous(labels=percent) + labs(y = "Percent of Total",x="Promoter width in base pairs (bp)")

ggsave("width_compare.png",width=6,height=5)

q()
