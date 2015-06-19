##creates a simple grid heatmap of meiosis genes
require(ggplot2)
require(reshape2)

meiosis_expression <- read.table(file="/home/rtraborn/Daphnia/Daphnia_CAGE_Data/R_analysis/expression_analysis/de_tables/meiosis_expression_table_de.txt",header=TRUE,stringsAsFactors=FALSE)
meiosis_table <- cbind(meiosis_expression$t_value_1,meiosis_expression$t_value_2,meiosis_expression$t_value_3) 
rownames(meiosis_table) <- rownames(meiosis_table)
colnames(meiosis_table) <- c("Ameiotic female - meiotic female","Meiotic female - male","Male - meiotic female")

plot.values <- melt(meiosis_table)
names(plot.values) <- c("gene","variable","value")    

p <- ggplot(plot.values, aes(variable, gene)) + geom_tile(aes(fill = value))

meiosis_figure <- p + scale_fill_gradient2(low="#D7191C", mid="white", high="#2C7BB6") +
labs(y="Gene symbol", x="", fill="t-value") +  scale_x_discrete(expand = c(0, 0)) +
     scale_y_discrete(expand = c(0, 0)) + theme(axis.text.x=element_text(angle = -45, hjust = 0))

png(file="meiosis_grid_heatmap.png",height=1200,width=900)
meiosis_figure


