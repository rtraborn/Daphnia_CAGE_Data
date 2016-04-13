
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

library(ggplot2)

library(reshape2)

setwd("/home/rtraborn/Daphnia/Daphnia_CAGE_Data/R_analysis/expression_analysis/")

de_results <- read.table(file="D_pulex_de_consensus_promoters.txt",header=TRUE)

head(de_results)

de_table <- array(NA,c(2,5))

colnames(de_table) <- colnames(de_results)

rownames(de_table) <- c("upregulated","downregulated")

head(de_table)

MvMF_up <- length(which(de_results$MvMf > 0))

MvMF_down <- length(which(de_results$MvMf < 0))

pEvM_up <- length(which(de_results$pEvM > 0))

pEvM_down <- length(which(de_results$pEvM < 0))

MfvpE_up <- length(which(de_results$MfvpE > 0))

MfvpE_down <- length(which(de_results$MfvpE < 0))

MvF_up <- length(which(de_results$MvF > 0))

MvF_down <- length(which(de_results$MvF < 0))

AvS_up <- length(which(de_results$AvS > 0))

AvS_down <- length(which(de_results$AvS < 0))

de_table[1,] <- c(MvMF_up,pEvM_up,MfvpE_up,MvF_up,AvS_up)

de_table[2,] <- c(-MvMF_down,-pEvM_down,-MfvpE_down,-MvF_down,-AvS_down)

de_table <- as.data.frame(de_table)

de_table$de <- c("upregulated","downregulated")

head(de_table)

de_table_m <- melt(de_table)

head(de_table_m)

p <- ggplot(subset(de_table_m,de=="upregulated"), aes(y=value,x=variable)) + geom_bar(stat="identity",colour="black",fill="yellow2") + theme_bw() + scale_y_continuous(limits = c(0,1250))

p

p2 <- ggplot(subset(de_table_m,de=="downregulated"), aes(y=value,x=variable)) + geom_bar(stat="identity",colour="black",fill="red3") + scale_y_continuous(limits = c(-1250,0))

p2 <- p2 + theme_bw() + theme(axis.title.x= element_blank(), axis.ticks.x = element_blank(),axis.text.x = element_blank())  

png(filename = "de_genes_by_comparison.png",width=7,height=7,units="in",res=300)

multiplot(p, p2, cols=1)

dev.off()

getwd()
