
library(gplots)

library(RColorBrewer)

setwd("/home/rtraborn/Daphnia/Daphnia_CAGE_Data/R_analysis/promoter_calling_pipelines/TCO/tagClusters/pooled_samples")

promoter_comp <- read.table(file="Dpm_core_matrix.logPvalue.matrix.txt", skip=1, header=FALSE,sep="\t",stringsAsFactors = FALSE)

colnames(promoter_comp) <- c("MotifID","Dpm1","Dpm2","Dpm3","Dpm4","Dpm5","Dpm6","Dpm7","Dpm8") 

row_in <- promoter_comp[,1]

promoter_comp <- promoter_comp[,-1]

rownames(promoter_comp) <- colnames(promoter_comp)

promoter_comp

promoter_comp_m <- as.matrix(promoter_comp)

head(promoter_comp_m)

is.matrix(promoter_comp_m)

r1 <- range(promoter_comp_m) - median(promoter_comp_m)

r1

hmcol<-brewer.pal(11,"RdBu")

par(mar=c(4.1,4.1,4.1,4.1))

png("Dpm_motifs_correlation.png",bg = "transparent",width= 1000, height = 1000, units = "px")

heatmap.2(promoter_comp_m,trace="none",notecol="black",col=colorRampPalette(c("red","white","blue"))(100))

dev.off()
