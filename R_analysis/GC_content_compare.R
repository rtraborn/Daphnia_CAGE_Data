library(ggplot2)
setwd("/home/rtraborn/Daphnia/Daphnia_CAGE_Data/R_analysis/promoter_calling_pipelines/TCO/tagClusters/pooled_samples")

Dp_random <- read.table(file="Dp_TCO_random_scaffold_intervals.nuc",skip=1,header=FALSE,stringsAsFactors=FALSE)
Dp_tag_clusters <- read.table(file="combined.tagClusters.qLow0.1_qUp0.9_merged.nuc",skip=1,header=FALSE,stringsAsFactors=FALSE)
Dp_tag_clusters[,7] -> TCs_vector
Dp_random[,8] -> random_vector
TC_len <- length(TCs_vector)
ran_len <- length(random_vector)
my_vec <- c(TCs_vector,random_vector)
combined_vector <- c(rep("Promoter",TC_len),rep("Background",ran_len))
TC_df <-cbind(my_vec,combined_vector)
TC_df2 <- data.frame(GC=as.numeric(TC_df[,1]),cond=factor(TC_df[,2]))

ggplot(TC_df2, aes(x=cond, y=GC, fill=cond)) + geom_boxplot() + guides(fill=FALSE)
ggsave(file="GC_content_compare.png",dpi=300,width=3,height=4)
