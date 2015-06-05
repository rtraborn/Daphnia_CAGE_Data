#creates a series of plots relating to differential expression analysis
# uses the bioconductor package 'limma'

require("CAGEr")
require("limma")
require("edgeR")
require("Biobase")
require("gplots")

#importing the Dp CAGEr object
load("/home/rtraborn/Daphnia/Daphnia_CAGE_Data/R_analysis/promoter_calling_pipelines/Dp_TCO.RData")

#CAGE data based expression clustering
getExpressionProfiles(myCAGEset, what = "consensusClusters",
                      tpmThreshold = 10, nrPassThreshold = 1,
                      method = "som", xDim = 4, yDim = 4
                      )


plotExpressionProfiles(myCAGEset, what= "consensusClusters")

#creating a count table
Dp_edger <- myCAGEset@consensusClustersTpmMatrix

#creating identifier for each tag cluster
Dp_edger_consensus_cluster <- consensusClusters(myCAGEset)
rownames(Dp_edger) <- paste(Dp_edger_consensus_cluster$chr,
					Dp_edger_consensus_cluster$start,
					Dp_edger_consensus_cluster$end,
					Dp_edger_consensus_cluster$strand,
					sep="_")



dp_eset <-new("ExpressionSet", exprs=as.matrix(Dp_edger))

#checking to see what the data.frame looks like
head(Dp_edger)

lib_sizes <- librarySizes(myCAGEset)
#contrast.matrix <- makeContrasts(mat_male-mat_fem, pE_fem-mat_male, pE_fem-mat_fem, levels=design)

#importing the design object
#load("Dp_design.RData")

design <- model.matrix(~ 0+factor(c(1,1,1,3,3,2,2,2)))
colnames(design) <- c("mat_fem", "pE_fem", "mat_male")
design

contrast.matrix <- makeContrasts(pE_fem-mat_fem, mat_male-pE_fem, mat_male-mat_fem, levels=design)

group <- c(rep("mat_fem",3),rep("mat_male",2),rep("pE_fem",3)) 
p_cutoff <- 0.01
rowsum_threshold <- 20

Dp_dge <- DGEList(counts=Dp_edger,group=group)
A <- rowSums(Dp_dge$counts) 
Dp_dge <- Dp_dge[A>rowsum_threshold,]  
Dp_dge <- calcNormFactors(Dp_dge)
Dp_dge <- estimateCommonDisp(Dp_dge, verbose=T)

de.tgw2_pE_male <- exactTest(Dp_dge,pair=c("pE_fem","mat_male"))
de.tgw2_pE_mat_fem <- exactTest(Dp_dge,pair=c("pE_fem","mat_fem"))
de.tgw2_male_mat_fem <- exactTest(Dp_dge,pair=c("mat_male","mat_fem"))

v <- voom(Dp_dge,design,plot=TRUE)
fit <- lmFit(v,design=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit)

topTable(fit2,coef=ncol(design))

volcanoplot(fit2)

options(digits=3)
top_table <- topTable(fit2,coef=ncol(design),n=Inf,sort.by="p",adjust="BH",p=0.01)
results <- decideTests(fit2)

vennDiagram(results)

#number of differentially-related promoters
sum(top_table$adj.P.Val<0.01)

de_data <- Dp_dge$pseudo.counts


#differential analysis results
de_data <- cbind(de_data, de.tgw2_pE_male$table)

#calculate FDR
de_data$FDR <- p.adjust(de_data$PValue, method='BH')

#dispersion of each tag cluster
de_data$tw_dis <- Dp_dge$tagwise.dispersion

#coordinates of each tag cluster
data_coord2 <- matrix(data=unlist(strsplit(rownames(de_data), split="_")),
                      nrow= length(row.names(de_data)),
                      byrow=T)
data_coord2 <- as.data.frame(data_coord2, stringsAsFactors=F)
head(data_coord2)

col_1 <- data_coord2[,1]
col_2 <- data_coord2[,2]
chr_col <- paste(col_1,col_2,sep="_")
data_coord2 <- data_coord2[,-2]
data_coord2[,1] <- chr_col
names(data_coord2) <- c('chr','start','end','strand')

#coordinates of each tag cluster
de_data <- cbind(de_data, data_coord2)

#create column for differential expression status
#1 for DE and 0 for not
de_data$de <- as.numeric(de_data$FDR<p_cutoff)

#convert coordinates to numeric
de_data$start <- as.numeric(de_data$start)
de_data$end <- as.numeric(de_data$end)

#creating a genomicRanges object
de_data_gr <- with(de_data, GRanges(chr,
                                    IRanges(start, end, names=row.names(de_data)), strand))

#what does the GR object look like?
de_data_gr

write.table(de_data,file="de_data_TCO.txt",col.names=TRUE,quote=FALSE)

###########################################################################

plotMA(fit2,values=c("mat_fem","pE_fem"),col=c("red","blue"),main="TCO MA Plot",legend="topright")

png(file="heatmap_TCO_all.png",height=900,width=1260)
selected  <- p.adjust(fit2$p.value[, 2]) <0.01
esetSel <- dp_eset[selected, ]
heatmap(exprs(esetSel))
dev.off()

#annotation (in progress- needs to be converted to Daphnia; this is a human example)
#columns <- c("GeneID","Symbol", "chromosome", "type_of_gene")
#GeneInfo <- read.columns("130811-Homo_sapiens.gene_info",required=columns, stringsAsFactors = FALSE)
#m <- match(gene$annotation$GeneID, GeneInfo$GeneID)

#Ann <- cbind(gene$annotation[, c("GeneID", "Chr", "Length")], GeneInfo[m, c("Symbol", "type_of_gene")])
#Ann$Chr <- unlist(lapply(strsplit(Ann$Chr, ";"), function(x) paste(unique(x), collapse = "|")))
#Ann$Chr <- gsub("chr", "", Ann$Chr)
#head(Ann)

