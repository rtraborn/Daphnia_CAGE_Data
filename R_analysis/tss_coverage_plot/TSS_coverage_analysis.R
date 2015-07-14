
setwd("/home/rtraborn/Daphnia/Daphnia_CAGE_Data/tss_coverage_plot/coverage_files")

library(ggplot2)

library(reshape2)

list.files(pattern = "*.coverage") -> coverage_files

options(scipen=999)

head(coverage_files)

my_annot <- list(length=8)

str(my_annot)

for (i in 1:length(coverage_files)){ read.table(file=coverage_files[i],header=FALSE) -> this_var; this_var -> my_annot[[i]] }

names(my_annot) <- coverage_files

head(str(my_annot))

count_vec <- vector(mode="list",length=length(coverage_files))

for (i in 1:length(coverage_files)) { my_annot[[i]] -> this_df; as.numeric(this_df$V7) -> this_count; sum(this_count) -> count_vec[i] }

my_annot[[8]] -> my_intergenic; sum(my_intergenic$V4) -> count_vec[8]

count_vec

cov_ids <- c("3_UTR_intron","3_UTR","500_upstream","5_UTR_intron","5_UTR","CDS_intron","CDS","Intergenic")

total_count <- sum(as.numeric(count_vec))

total_count

fraction_vec <- vector(mode="list",length=8)

for (i in 1:8) { as.numeric(count_vec[i])/total_count -> this_fraction; this_fraction -> fraction_vec[i]}

fraction_vec

fraction_df <- data.frame(cov_ids,as.numeric(fraction_vec),stringsAsFactors=FALSE)

is(fraction_df)

fraction_df$Sample <- "total"

fraction_df$Sample <- as.factor(fraction_df$Sample)

fraction_df[,1] <- factor(c("3_UTR_intron","3_UTR","500_upstream","5_UTR_intron","5_UTR","CDS_intron","CDS","Intergenic"),levels=c("500_upstream","5_UTR","CDS","Intergenic","CDS_intron","5_UTR_intron","3_UTR","3_UTR_intron"))

fraction_df

colnames(fraction_df) <- c("genome_region","coverage_fraction","Sample")

str(fraction_df)

fraction_df <- fraction_df[order(as.numeric(fraction_df$coverage_fraction),decreasing=TRUE) ,]

fraction_df

my_cov_plot <- ggplot(fraction_df, aes(y=coverage_fraction,x=Sample,fill=genome_region)) + geom_bar(width=0.40,height=10,stat = "identity")

my_cov_plot + coord_flip() + theme(legend.position="none")

ggsave("tss_coverage_plot.png",dpi=300,width=7,height=2)
