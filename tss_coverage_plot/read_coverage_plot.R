#written on April 8, 2015

library("ggplot2")

#importing the coverage data
total_cov <- read.table(file="total_read_coverage_i.txt",header=TRUE)
total_cov <- as.data.frame(total_cov)
total_cov$sample <- "total"
names(total_cov) <- c("coverage_fraction","region","sample")

#setting the colors (4)
mycols <- c('#FFFD00', '#97CB00', '#3168FF', '#FF0200')

#plotting the figure
png("TSS_coverage_plot.png",width=600,height=100)
ggplot(total_cov, aes(y=coverage_fraction,x=sample,fill=region))+geom_bar(width=0.40,stat='identity')+coord_flip()
dev.off()

