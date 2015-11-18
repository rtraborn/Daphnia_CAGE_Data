##creates a simple grid heatmap of meiosis genes
require(ggplot2)
require(reshape2)

meiosis_expression <- read.csv(file="/home/rtraborn/Daphnia/Daphnia_CAGE_Data/R_analysis/expression_analysis/de_tables/meiosis_expression_table_de.csv",header=TRUE,stringsAsFactors=FALSE)
meiosis_table <- cbind(meiosis_expression$t_value_Male_v_Females,meiosis_expression$t_value_Asex_v_Sexuals,meiosis_expression$t_value_Asex_v_pE) 
rownames(meiosis_table) <- meiosis_expression[,1]
colnames(meiosis_table) <- c("Male_v_Females","Asex_v_sexuals","Asex_v_pE")

plot.values <- melt(meiosis_table)
names(plot.values) <- c("gene","variable","value")    

p <- ggplot(plot.values, aes(variable, gene)) + geom_tile(aes(fill = value))

meiosis_figure <- p + scale_fill_gradient2(low="#D7191C", mid="white", high="#2C7BB6") +
labs(y="Gene symbol", x="", fill="t-value") +  scale_x_discrete(expand = c(0, 0)) +
     scale_y_discrete(expand = c(0, 0)) + theme(axis.text.x=element_text(angle = -45, hjust = 0))

png(file="meiosis_grid_heatmap.png",height=1200,width=900)
meiosis_figure


