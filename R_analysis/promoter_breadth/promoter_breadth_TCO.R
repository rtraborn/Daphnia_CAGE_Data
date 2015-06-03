###calculates promoter breadth in Daphnia and Drosophila using CAGE annotations

require("GenomicRanges")
require("ggplot2")
require("scales")

dpul.df <- read.table(file="/home/rtraborn/Daphnia/Daphnia_CAGE_Data/data_tables/TCO/TCO_combined_cons_table.txt",header=TRUE,stringsAsFactors=FALSE)

names(dpul.df) <- c("chr","start","end","strand","geneID","start_gene","end_gene","strand2","tpm")
                                                                                                                              dpul.gr <- makeGRangesFromDataFrame(dpul.df, seqnames.field="chr", start.field="start", end.field="end", strand.field="strand")

dpul_widths <- width(dpul.gr)
n_widths_dp <- length(dpul_widths)

width_vec <- c(dpul_widths)
org_vec <- c(rep("D_pulex",n_widths_dp))

widths.df <- data.frame(width_vec, org_vec)
names(widths.df) <- c("promoter_breadth","species")

#plotting the figures

ggplot(widths.df, aes(x=promoter_breadth,y = (..count..)/sum(..count..))) + geom_histogram(binwidth=5,colour="orange2",fill="white") + facet_grid(species ~ .) + scale_y_continuous(labels=percent) + labs(y = "Percent of Total",x="Promoter width in base pairs (bp)")

ggsave("width_TCO.png",width=8,height=6)

q()
