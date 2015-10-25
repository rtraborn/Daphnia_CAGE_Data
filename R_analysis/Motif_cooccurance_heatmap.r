
library(gplots)

library(RColorBrewer)

setwd("/home/rtraborn/Daphnia/Daphnia_CAGE_Data/R_analysis/promoter_calling_pipelines/TCO/tagClusters/pooled_samples/homer_pos_files")

promoter_comp <- read.table(file="DpTCO_occ.logPvalue.matrix.txt", skip=1, header=FALSE,sep="\t",stringsAsFactors = FALSE)

head(promoter_comp)

colnames(promoter_comp) <- c("MotifID","Dpm1","Dpm2","Dpm3","Dpm4","Dpm5","Dpm6","Dpm7") #renaming to match the final Dpm set

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

png(file="Dpm_motifs_correlation.png",height=960,width=960,bg = "transparent")

plot.new()

heatmap.2(promoter_comp_m,trace="none",notecol="black",col=colorRampPalette(c("red","white","blue"))(100))

dev.off()

list.files()

DpTCO_peak_Dpm <- read.table(file="DpTCO_peak_Dpm_core_TCs_combined_sorted.txt",header=TRUE,sep="\t",na.strings=c("", "NA"))

head(DpTCO_peak_Dpm)

DpTCO_breadth_SI <- read.table(file="~/Daphnia//Daphnia_CAGE_Data/R_analysis/promoter_calling_pipelines/TCO/tagClusters/pooled_samples/TSR_breadth_SI_9_2_15.txt",header=TRUE,sep="\t")

head(DpTCO_breadth_SI)

cbind(DpTCO_breadth_SI,DpTCO_peak_Dpm[,6:29]) -> DpTCO_motif_breadth_combined

head(DpTCO_motif_breadth_combined)

write.table(DpTCO_motif_breadth_combined,file="DpTCO_motif_breadth_combined.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote = FALSE)

mean(as.numeric(na.omit(DpTCO_motif_breadth_combined$CpG.)))

mean(as.numeric(na.omit(DpTCO_motif_breadth_combined$GC.)))

Dpm1_occur <- DpTCO_motif_breadth_combined$X...Dpm1.Distance.From.Peak.sequence.strand.conservation.

length(na.omit(Dpm1_occur))

Dpm2_occur <- DpTCO_motif_breadth_combined$Dpm2.Distance.From.Peak.sequence.strand.conservation.

length(na.omit(Dpm2_occur))

Dpm3_occur <- DpTCO_motif_breadth_combined$Dpm3.Distance.From.Peak.sequence.strand.conservation.

length(na.omit(Dpm3_occur))

Dpm4_occur <- DpTCO_motif_breadth_combined$Dpm4.Distance.From.Peak.sequence.strand.conservation.

length(na.omit(Dpm4_occur))

Dpm5_occur <- DpTCO_motif_breadth_combined$Dpm5.Distance.From.Peak.sequence.strand.conservation.

length(na.omit(Dpm5_occur))

Dpm6_occur <- DpTCO_motif_breadth_combined$Dpm6.Distance.From.Peak.sequence.strand.conservation.

length(na.omit(Dpm6_occur))

Dpm7_occur <- DpTCO_motif_breadth_combined$Dpm7.Distance.From.Peak.sequence.strand.conservation.

length(na.omit(Dpm7_occur))

Dpm8_occur <- DpTCO_motif_breadth_combined$Dpm8.Distance.From.Peak.sequence.strand.conservation.

Dpm_NA_index <- which(is.na(DpTCO_motif_breadth_combined$X...Dpm1.Distance.From.Peak.sequence.strand.conservation.))

Dpm_presence_array <- array(NA,c(nrow(DpTCO_motif_breadth_combined),8)) 

c(0) -> Dpm_presence_array[Dpm_NA_index,1]; c(1) -> Dpm_presence_array[-Dpm_NA_index,1]

Dpm_NA_index_2 <- which(is.na(DpTCO_motif_breadth_combined$Dpm2.Distance.From.Peak.sequence.strand.conservation.))

c(0) -> Dpm_presence_array[Dpm_NA_index_2,2]; c(1) -> Dpm_presence_array[-Dpm_NA_index_2,2]

Dpm_NA_index_3 <- which(is.na(DpTCO_motif_breadth_combined$Dpm3.Distance.From.Peak.sequence.strand.conservation.))

c(0) -> Dpm_presence_array[Dpm_NA_index_3,3]; c(1) -> Dpm_presence_array[-Dpm_NA_index_3,3]

Dpm_NA_index_4 <- which(is.na(DpTCO_motif_breadth_combined$Dpm4.Distance.From.Peak.sequence.strand.conservation.))

c(0) -> Dpm_presence_array[Dpm_NA_index_4,4]; c(1) -> Dpm_presence_array[-Dpm_NA_index_4,4]

Dpm_NA_index_5 <- which(is.na(DpTCO_motif_breadth_combined$Dpm5.Distance.From.Peak.sequence.strand.conservation.))

c(0) -> Dpm_presence_array[Dpm_NA_index_5,5]; c(1) -> Dpm_presence_array[-Dpm_NA_index_5,5]

Dpm_NA_index_6 <- which(is.na(DpTCO_motif_breadth_combined$Dpm6.Distance.From.Peak.sequence.strand.conservation.))

c(0) -> Dpm_presence_array[Dpm_NA_index_6,6]; c(1) -> Dpm_presence_array[-Dpm_NA_index_6,6]

Dpm_NA_index_7 <- which(is.na(DpTCO_motif_breadth_combined$Dpm7.Distance.From.Peak.sequence.strand.conservation.))

c(0) -> Dpm_presence_array[Dpm_NA_index_7,7]; c(1) -> Dpm_presence_array[-Dpm_NA_index_7,7]

Dpm_NA_index_8 <- which(is.na(DpTCO_motif_breadth_combined$Dpm8.Distance.From.Peak.sequence.strand.conservation.))

c(0) -> Dpm_presence_array[Dpm_NA_index_8,8]; c(1) -> Dpm_presence_array[-Dpm_NA_index_8,8]

colnames(Dpm_presence_array) <- (c("Dpm1","Dpm2","Dpm3","Dpm4","Dpm5","Dpm6","Dpm7","Dpm8"))

head(Dpm_presence_array)

cbind(DpTCO_breadth_SI,Dpm_presence_array) -> Dpm_breadth_motif_presence_matrix

head(Dpm_breadth_motif_presence_matrix)

write.table(Dpm_breadth_motif_presence_matrix,file="Dpm_breadth_motif_presence_matrix.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote = FALSE)

summary(Dpm_breadth_motif_presence_matrix)

peaked_index <- which(Dpm_breadth_motif_presence_matrix$SI>=0)

broad_index <- which(Dpm_breadth_motif_presence_matrix$SI<=-1)

shape <- rep(NA,nrow(Dpm_breadth_motif_presence_matrix))

shape[peaked_index] <- "peaked"

shape[broad_index] <- "broad"

combined_index <- c(peaked_index,broad_index)

head(combined_index)

shape[-combined_index] <- "unclassified"

Dpm_breadth_motif_presence_matrix_i <- cbind(Dpm_breadth_motif_presence_matrix,shape)

promoter_label <- seq(1:nrow(Dpm_breadth_motif_presence_matrix_i))

promoter_ID <- paste("promoter",promoter_label,sep="_")

Dpm_breadth_motif_presence_matrix_i <- cbind(promoter_ID, Dpm_breadth_motif_presence_matrix_i)

head(Dpm_breadth_motif_presence_matrix_i)

library(reshape2)

library(ggplot2)

df_1 <- cbind(as.character(Dpm_breadth_motif_presence_matrix_i$promoter_ID), Dpm_breadth_motif_presence_matrix_i$breadth, Dpm_breadth_motif_presence_matrix_i$nTSSs, Dpm_breadth_motif_presence_matrix_i$SI, as.character(Dpm_breadth_motif_presence_matrix_i$shape),Dpm_breadth_motif_presence_matrix_i$Dpm2,Dpm_breadth_motif_presence_matrix_i$Dpm3)

head(df_1)

df_1 <- as.data.frame(df_1)

colnames(df_1) <- c("promoter_name","breadth","nTSSs","SI","Shape","TATA","Initiator")

head(df_1)

df_melt <- melt(df_1,varnames=c(Shape,SI,nTSSs,TATA,Initiator))

head(df_melt)

df_melt <- df_melt[,-1]

ggplot(df_melt, aes(x=factor(Shape),y=as.numeric(nTSSs))) + geom_boxplot()

df_melt2 <- df_melt[,-4]

head(df_melt2)

TATA_index <- which(df_melt2[,4]==1)

TATA_df <- df_melt[TATA_index,]

Inr_index <- which(df_melt2[,5]==1)

Inr_df <- df_melt[Inr_index,]

head(Inr_df)

ggplot(Inr_df, aes(x=factor(Shape),y=as.numeric(breadth))) + geom_boxplot()

ggplot(Inr_df, aes(x=factor(Shape),y=as.numeric(nTSSs))) + geom_boxplot()

ggplot(TATA_df, aes(x=factor(Shape),y=as.numeric(breadth))) + geom_boxplot()

combined_df <- rbind(TATA_df,Inr_df) 

nrow(Inr_df)

motif_vector <- c(rep("TATA",900),rep("Inr",3392))

combined_df <- cbind(combined_df,motif_vector)

head(combined_df)

#ggplot(combined_df, aes(x=factor(motif_vector),y=as.numeric(SI))) + geom_boxplot()

data_array <- array(NA, c(nrow(df_melt2),5))

data_array[,1] <- as.numeric(df_melt2[,1])

data_array[,2] <- as.numeric(df_melt2[,2])

data_array[,3] <- as.numeric(df_melt2[,3])

data_array[,4] <- as.numeric(df_melt2[,4])

data_array[,5] <- as.numeric(df_melt2[,5])

colnames(data_array) <- colnames(df_melt2)

cormat <- cor(data_array)

head(cormat)

sum(data_array[,4])

sum(data_array[,5])

nrow(data_array)

head(Dpm_breadth_motif_presence_matrix_i)

Dpm_presence_array_i <- Dpm_presence_array[,-6]

head(Dpm_presence_array_i)

sum(Dpm_presence_array_i[,1])

my_matrix <- data.matrix(Dpm_presence_array_i)

head(my_matrix)

png(file="motif_cooc_heatmap_all.png")

#heatmap.2(Dpm_presence_array_i,scale="none")

TATA_index <- which(Dpm_breadth_motif_presence_matrix_i$Dpm2==1)

head(data_array)

head(combined_df)

dim(combined_df)

dim(Dpm_breadth_motif_presence_matrix_i)

head(Dpm_breadth_motif_presence_matrix_i)

TATA <- rep("TATA-less",nrow(Dpm_breadth_motif_presence_matrix_i))

TATA[TATA_index] <- "TATA"

head(TATA_vector)

length(TATA_index)

Dpm_total_matrix <- cbind(Dpm_breadth_motif_presence_matrix_i,TATA)

head(Dpm_total_matrix)

ggplot(Dpm_total_matrix, aes(x=factor(TATA),y=as.numeric(SI))) + geom_boxplot()

head(Dpm_breadth_motif_presence_matrix_i)

new_matrix <- cbind(Dpm_breadth_motif_presence_matrix_i$breadth,Dpm_breadth_motif_presence_matrix_i$nTSSs, Dpm_breadth_motif_presence_matrix_i$SI)

colnames(new_matrix) <- c("breadth","count","SI")

new_df <- as.data.frame(new_matrix)

head(new_df)

dim(new_matrix)

quantile(new_df$breadth,probs = seq(0, 1, 0.05))

quantile(new_df$SI,probs = seq(0, 1, 0.05))

broadest_index <- which(new_df$breadth > 45)

sharpest_index <- which(new_df$SI > 1.12382271532844)

length(sharpest_index)

breadth_id <- rep("regular",nrow(new_matrix))

breadth_id[broadest_index] <- "broad"

breadth_id[sharpest_index ] <- "peaked"

new_df2 <- cbind(new_df,breadth_id)

head(new_df2)

ggplot(new_df2, aes(x=count, y=breadth, color=breadth_id)) + geom_point(shape=1) + scale_x_continuous(limit=c(0,10000)) + geom_smooth(method=lm) 

test_loess <- loess(data=new_df2,count ~ breadth)

summary(test_loess)

test_lm <- lm(data=new_df2,SI ~ breadth)

summary(test_lm)

ggplot(new_df2, aes(x=count, y=breadth, colour=breadth_id)) + geom_point(shape=1) + scale_x_log10() + geom_smooth(model="loess")

ggsave(file="consensus_TSR_shape_xy_scatter.png")

head(new_df2)

test_lm <- lm(formula = count ~ breadth, data=new_df2)

summary(test_lm)

ggplot(new_df2, aes(x=breadth_id, y=count,colour=breadth_id)) + geom_boxplot() + scale_y_continuous(limit=c(0,50000))

ggsave(file="consensus_TSR_shape_boxplot.png")
