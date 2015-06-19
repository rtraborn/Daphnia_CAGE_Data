#written on April 8, 2015

library("ggplot2")
library("reshape2")

#importing the coverage data
total_cov <- read.table(file="total_read_coverage.txt",header=TRUE)

#plotting the figure
png("TSS_coverage_plot.png",width=960,height=960)
ggplot(total_cov, aes(x = stage, y = fraction_of_reads,fill=region))+geom_bar(width=0.75,stat='identity')+coord_flip()
dev.off()
