
library(gplots)

library(RColorBrewer)

setwd("/home/rtraborn/Daphnia/Daphnia_CAGE_Data/R_analysis/promoter_calling_pipelines/TCO/tagClusters/pooled_samples/homer_pos_files")

promoter_comp <- read.table(file="DpTCO_occ.logPvalue.matrix.txt", skip=1, header=FALSE,sep="\t",stringsAsFactors = FALSE)

head(promoter_comp)

colnames(promoter_comp) <- c("MotifID","Motif1","Motif2","Motif3","Motif4","Motif5","Motif6","Motif7","Motif8") #renaming to match the final Dpm set

row_in <- promoter_comp[,1]

promoter_comp <- promoter_comp[,-1]

head(promoter_comp)

rownames(promoter_comp) <- colnames(promoter_comp)

promoter_comp

promoter_comp_m <- as.matrix(promoter_comp)

head(promoter_comp_m)

is.matrix(promoter_comp_m)

r1 <- range(promoter_comp_m) - median(promoter_comp_m)

r1

hmcol<-brewer.pal(11,"RdBu")

png(file="Dpm_motifs_co_occurance.png",height=960,width=960,bg = "transparent")

plot.new()

heatmap.2(promoter_comp_m,trace="none",notecol="black",col=colorRampPalette(c("red","white","blue"))(100))

dev.off()
