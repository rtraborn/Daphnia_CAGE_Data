#creates a series of plots relating to differential expression analysis
# uses the bioconductor package 'limma'

require("CAGEr")
require("limma")
require("edgeR")
require("Biobase")
require("gplots")
require(GenomicRanges)
require(GenomicFeatures)

#importing the Dp CAGEr object
load("/home/rtraborn/Daphnia/Daphnia_CAGE_Data/R_analysis/promoter_calling_pipelines/TCO/Dp_CAGEset_TCO.RData")

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

head(dp_eset)

#checking to see what the data frame looks like
head(Dp_edger)

lib_sizes <- librarySizes(myCAGEset)
design <- model.matrix(~0+factor(c(1,1,1,3,3,2,2,2)))
colnames(design) <- c("mat_fem", "pE_fem", "mat_male")
design

contrasts.matrix <- makeContrasts(
    MvMf=mat_male-mat_fem,
    pEvM=pE_fem-mat_male,
    MfvpE=mat_fem-pE_fem,
    MvF=mat_male-(mat_fem+pE_fem)/2,
    AvS=mat_fem-(mat_male+pE_fem)/2,
    levels=design)

#for direct comparison between the three groups
#contrasts.matrix <- makeContrasts(pE_fem-mat_fem, mat_male-pE_fem, mat_male-mat_fem,
#                                  levels=design)

contrasts.matrix

group <- c(rep("mat_fem",3),rep("mat_male",2),rep("pE_fem",3))

p_cutoff <- 0.01
rowsum_threshold <- 25

Dp_dge <- DGEList(counts=Dp_edger,group=group)
A <- rowSums(Dp_dge$counts) 
Dp_dge <- Dp_dge[A>rowsum_threshold,]  
Dp_dge <- calcNormFactors(Dp_dge)
Dp_dge <- estimateCommonDisp(Dp_dge, verbose=T)
Dp_dge <- estimateTagwiseDisp(Dp_dge, trend="none")

plotBCV(Dp_dge)

#de.tagwisedisp <- exactTest(Dp_dge)

v <- voom(Dp_dge, design, plot=TRUE)
fit <- lmFit(v,design)
fit2 <- contrasts.fit(fit,contrasts.matrix)
fit2 <- eBayes(fit2)
results <- decideTests(fit2,p.value=0.01,lfc=log2(2))
#ind = which(apply(res,1,function(x) {length(which(x != 0))>0}) == T)
#length(ind)

head(res)

options(digits=3)
de_table1 <- topTable(fit2, coef=1, sort="none",adjust.method="BH",n=Inf)
head(de_table1)
de_table2 <- topTable(fit2, coef=2, sort="none",adjust.method="BH",n=Inf)
head(de_table2)
de_table3 <- topTable(fit2, coef=3, sort="none",adjust.method="BH",n=Inf)
head(de_table3)
de_table4 <- topTable(fit2, coef=4, sort="none", adjust.method="BH",n=Inf)
de_table5 <- topTable(fit2, coef=5, sort="none", adjust.method="BH",n=Inf)

plotMDS(v, labels=c("mf1","mf2","mf3","m1","m2","pE1","pE2","pE3"), main="MDS plot for all eight libraries")
volcanoplot(fit2,coef=1,highlight=20)
volcanoplot(fit2,coef=2,highlight=20)
volcanoplot(fit2,coef=3,highlight=20)
volcanoplot(fit2,coef=4,highlight=20)
volcanoplot(fit2,coef=5,highlight=20)

#vennDiagram(results, names=c("Asexual females","Sexual females","Sexual Males"),include="both",circle.col=c("green","red","blue"))
#vennDiagram(results, names=c("Asexual females","Sexual females","Sexual Males"),include="up",circle.col=c("green","red","blue"))
#vennDiagram(results, names=c("Asexual females","Sexual females","Sexual Males"),include="down",circle.col=c("green","red","blue"))

top_table2 <- cbind(top_table,results)
write.table(top_table2,file="de_TCO_combined_table.txt", row.names=TRUE,quote=FALSE)

################################### mat females vs pE females ############################

de_data <- Dp_dge$pseudo.counts

#differential analysis results
de_data <- cbind(de_data, de_table1)

#calculating the false discovery rate
de_data$FDR <- p.adjust(de_data$P.Value, method = 'BH')

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

de_data1 <- de_data

write.table(de_data1,file="TCO_male_v_mat_fem_de.txt",col.names=TRUE,quote=FALSE)

################################### pE females vs males ############################

de_data <- Dp_dge$pseudo.counts

#differential analysis results
de_data <- cbind(de_data, de_table2)

de_data$FDR <- p.adjust(de_data$P.Value, method = 'BH')

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

de_data2 <- de_data

write.table(de_data2,file="TCO_pE_v_males_de.txt",col.names=TRUE,quote=FALSE)

###################################  asexual females vs pE_fem  ############################

de_data <- Dp_dge$pseudo.counts

#differential analysis results
de_data <- cbind(de_data, de_table3)

#calculating the FDR
de_data$FDR <- p.adjust(de_data$P.Value, method = 'BH')

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

de_data3 <- de_data

write.table(de_data3,file="TCO_matfem_v_pE_de.txt",col.names=TRUE,quote=FALSE)
###################################  Males vs Females  ############################

de_data <- Dp_dge$pseudo.counts

#differential analysis results
de_data <- cbind(de_data, de_table4)

#calculating the FDR
de_data$FDR <- p.adjust(de_data$P.Value, method = 'BH')

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

de_data4 <- de_data

write.table(de_data4,file="TCO_male_v_females_de.txt",col.names=TRUE,quote=FALSE)

###################################  Asexual vs Sexuals  ############################

de_data <- Dp_dge$pseudo.counts

#differential analysis results
de_data <- cbind(de_data, de_table5)

#calculating the FDR
de_data$FDR <- p.adjust(de_data$P.Value, method = 'BH')

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

de_data5 <- de_data

write.table(de_data5,file="TCO_asexual_v_sexuals_de.txt",col.names=TRUE,quote=FALSE)

####################################

#Adding gene annotation to promoters

#importing the gene annotation file (only need to do this a single time)
dpulex_genes <- read.table(file="/home/rtraborn/Daphnia/Daphnia_CAGE_Data/gene_annotations/dpulex_genes.bed", header=FALSE)
names(dpulex_genes) <- c("chr","start","end","geneID","score","strand","version","type","placeholder","ID2")

#create GRanges object for all genes
genes_GR <- with(dpulex_genes, GRanges(chr, IRanges(start, end, names=geneID),strand))

#################################### pE vs males ########################################

#next, creating a genomicRanges object from the promoter data
de_GR <- with(de_data1, GRanges(chr,
                                    IRanges(start, end, names=row.names(de_data1)), strand))

de_GR <- promoters(de_GR, upstream=500, downstream=100)

#de_GR #for debugging

#overlap TSRs with all tag clusters
tc_overlaps <- findOverlaps(de_GR, genes_GR)

#store number of overlaps
#tc_overlaps_count <- countOverlaps(de_GR, genes_GR)

match_hit2 <- as.data.frame(tc_overlaps)

#name the columns
names(match_hit2) <- c('query','subject')

#head(match_hit2)

promoter_index <- match_hit2$query
gene_index <- match_hit2$subject

print("Length of promoter index")
length(promoter_index)

#gene_names <- dpulex_genes[gene_index,'geneID']
promoter_table <- data.frame(de_data1[promoter_index,])
gene_IDs <- dpulex_genes[gene_index,"geneID"]
promoter_table$gene <- gene_IDs

#remove duplicated entries
#match_hit2 <- match_hit2[!duplicated(match_hit2$query),]

de_index <- which(promoter_table$FDR<0.01)
length(de_index)
promoter_table <- promoter_table[de_index,]

write.table(promoter_table,file="TCO_MvMf_de_table_genes.txt",col.names=TRUE, row.names=TRUE)

#################################### pE vs asexual females  ########################################

#creating a genomicRanges object from the promoter data
de_GR <- with(de_data2, GRanges(chr,
                                    IRanges(start, end, names=row.names(de_data2)), strand))

de_GR <- promoters(de_GR, upstream=200, downstream=200)

#de_GR #for debugging

#overlap TSRs with all tag clusters
tc_overlaps <- findOverlaps(de_GR, genes_GR)

#store number of overlaps
#tc_overlaps_count <- countOverlaps(de_GR, genes_GR)

match_hit2 <- as.data.frame(tc_overlaps)

#name the columns
names(match_hit2) <- c('query','subject')

#head(match_hit2)

promoter_index <- match_hit2$query
gene_index <- match_hit2$subject

print("Length of promoter index")
length(promoter_index)

#gene_names <- dpulex_genes[gene_index,'geneID']
promoter_table <- data.frame(de_data2[promoter_index,])
gene_IDs <- dpulex_genes[gene_index,"geneID"]
promoter_table$gene <- gene_IDs

#remove duplicated entries
#match_hit2 <- match_hit2[!duplicated(match_hit2$query),]

de_index <- which(promoter_table$FDR<0.01)
length(de_index)
promoter_table <- promoter_table[de_index,]

write.table(promoter_table,file="TCO_pEvM_de_table_genes.txt",col.names=TRUE, row.names=TRUE)

#################################### pE vs asexual females  ########################################

#creating a genomicRanges object from the promoter data
de_GR <- with(de_data3, GRanges(chr,
                                    IRanges(start, end, names=row.names(de_data3)), strand))

de_GR <- promoters(de_GR, upstream=200, downstream=200)

#de_GR #for debugging

#overlap TSRs with all tag clusters
tc_overlaps <- findOverlaps(de_GR, genes_GR)

#store number of overlaps
#tc_overlaps_count <- countOverlaps(de_GR, genes_GR)

match_hit2 <- as.data.frame(tc_overlaps)

#name the columns
names(match_hit2) <- c('query','subject')

#head(match_hit2)

promoter_index <- match_hit2$query
gene_index <- match_hit2$subject

print("Length of promoter index")
length(promoter_index)

#gene_names <- dpulex_genes[gene_index,'geneID']
promoter_table <- data.frame(de_data3[promoter_index,])
gene_IDs <- dpulex_genes[gene_index,"geneID"]
promoter_table$gene <- gene_IDs

#remove duplicated entries
#match_hit2 <- match_hit2[!duplicated(match_hit2$query),]

de_index <- which(promoter_table$FDR <0.01)
length(de_index)
promoter_table <- promoter_table[de_index,]

write.table(promoter_table,file="TCO_MfvpE_de_table_genes.txt",col.names=TRUE, row.names=TRUE)

####################################

#number of differentially-related promoters
sum(de_data1$FDR<0.01)
sum(de_data2$FDR<0.01)
sum(de_data3$FDR<0.01)
#sum(de_table$B<0.01)

################################# Gene Family Analyses ##############################

#importing the meiosis gene annotation file
meiosis_genes <- read.table(file="/home/rtraborn/Daphnia/Daphnia_CAGE_Data/development_reg/meiosis/meiosis_genes_upstream_intervals.bed", header=TRUE)
meiosis_IDs <- meiosis_genes$geneID
promoter_list <- promoter_table$gene

meiosis_gr <- with(meiosis_genes, GRanges(chr,
                                    IRanges(start, end, names=meiosis_IDs), strand))

meiosis_OL <- findOverlaps(de_GR, meiosis_gr)

match_hit2 <- as.data.frame(meiosis_OL)

#name the columns
names(match_hit2) <- c('query','subject')

#head(match_hit2)

promoter_index <- match_hit2$query
gene_index <- match_hit2$subject

print("Length of promoter index")
length(promoter_index)

#gene_names <- dpulex_genes[gene_index,'geneID']
promoter_table <- data.frame(de_data1[promoter_index,])
gene_IDs <- dpulex_genes[gene_index,"geneID"]
promoter_table$gene <- gene_IDs

#remove duplicated entries
#match_hit2 <- match_hit2[!duplicated(match_hit2$query),]

meiosis_table <- na.omit(promoter_table)

#differentially-expressed meiosis genes only
de_index <- which(meiosis_table$de == 1)
meiosis_de <- meiosis_table[de_index,]

#head(meiosis_table)
write.table(meiosis_table,file="meiosis_table.txt",col.names=TRUE,row.names=TRUE,quote=FALSE)
write.table(meiosis_de,file="meiosis_table_de.txt",col.names=TRUE,row.names=TRUE,quote=FALSE)

##########################################################################
#importing the gust. receptors  gene annotation file
gust_genes <- read.table(file="/home/rtraborn/Daphnia/Daphnia_CAGE_Data/development_reg/gustatory_rec/Dpulex_gustatory_receptor_genes.bed", header=TRUE,stringsAsFactors=FALSE)
rownames(gust_genes) <- gust_genes$geneID
gust_IDs <- gust_genes$geneID
head(gust_IDs)
promoter_list <- unlist(promoter_table$gene)
promoter_list <- as.character(promoter_list)
head(promoter_list)
my_index <-  match(gust_IDs,promoter_list)
length(my_index)
head(my_index)
gustatory_table <- promoter_table[my_index,]
gustatory_table <- na.omit(gustatory_table)
head(gustatory_table)

#differentially-expressed gustatory genes only
de_index <- which(gustatory_table$de == 1)
gustatory_de <- gustatory_table[de_index,]

#write.table(gustatory_table,file="gustatory_table.txt",col.names=TRUE,row.names=TRUE,quote=FALSE)
#write.table(gustatory_de,file="gustatory_table_de.txt",col.names=TRUE,row.names=TRUE,quote=FALSE)

###########################################################################
#Making heatmaps from the eset data we've generated

#overall heatmap
par(mar=c(4.1,4.1,4.1,4.1))
png(file="all_genes_heatmap.png",height=2800,width=2800)
selected  <- rownames(de_data1)
esetSel <- dp_eset[selected,]
heatmap.2(exprs(esetSel), symm=FALSE,symkey=FALSE,scale="row", density.info="none",trace="none",
          key=TRUE,margins=c(10,10))
dev.off()

#male vs asexuals
par(mar=c(4.1,4.1,4.1,4.1))
png(file="heatmap_TCO_male_v_matfem.png",height=2800,width=2800)
de_index <- which(de_data1$FDR<0.01)
length(de_index)
de <- de_table1[de_index,]
selected <- rownames(de)
esetSel <- dp_eset[selected, ]
heatmap.2(exprs(esetSel), symm=FALSE,symkey=FALSE,scale="row", density.info="none",trace="none",
          key=TRUE,margins=c(10,10))
dev.off()

#pE v male
par(mar=c(4.1,4.1,4.1,4.1))
png(file="heatmap_TCO_pE_v_male.png",height=2800,width=2800)
de_index <- which(de_data2$FDR<0.01)
length(de_index)
de <- de_table2[de_index,]
selected <- rownames(de)
esetSel <- dp_eset[selected, ]
heatmap.2(exprs(esetSel), symm=FALSE,symkey=FALSE,scale="row", density.info="none",trace="none",
          key=TRUE,margins=c(10,10))
dev.off()

#asexuals vs pE
par(mar=c(4.1,4.1,4.1,4.1))
png(file="heatmap_TCO_matfem_v_pE.png",height=2800,width=2800)
de_index <- which(de_data3$FDR<0.01)
length(de_index)
de <- de_table3[de_index,]
selected <- rownames(de)
esetSel <- dp_eset[selected, ]
heatmap.2(exprs(esetSel), symm=FALSE,symkey=FALSE, scale="row", density.info="none",trace="none",
         key=TRUE, margins=c(10,10))
dev.off()

#male vs both females
par(mar=c(4.1,4.1,4.1,4.1))
png(file="heatmap_TCO_males_v_females.png",height=2800,width=2800)
de_index <- which(de_data4$FDR<0.01)
length(de_index)
de <- de_table4[de_index,]
selected <- rownames(de)
esetSel <- dp_eset[selected, ]
heatmap.2(exprs(esetSel), symm=FALSE,symkey=FALSE, scale="row", density.info="none",trace="none",
         key=TRUE, margins=c(10,10))
dev.off()

#asexuals vs both sexuals
par(mar=c(4.1,4.1,4.1,4.1))
png(file="heatmap_asex_v_sexuals.png",height=2800,width=2800)
de_index <- which(de_data5$FDR<0.01)
length(de_index)
de <- de_table5[de_index,]
selected <- rownames(de)
esetSel <- dp_eset[selected, ]
heatmap.2(exprs(esetSel), symm=FALSE,symkey=FALSE, scale="row", density.info="none",trace="none",
         key=TRUE, margins=c(10,10))
dev.off()

#meiosis genes
par(mar=c(4.1,4.1,4.1,4.1))
png(file="heatmap_TCO_meiosis.png",height=2800,width=2800)
meiosis_rows <- match(rownames(meiosis_table), rownames(de_table1))
length(meiosis_rows)
head(meiosis_rows)
meiosis_rows <- na.omit(meiosis_rows)
selected  <- rownames(de_data1[meiosis_rows,])
esetSel <- dp_eset[selected, ]
heatmap.2(exprs(esetSel), symm=FALSE,symkey=FALSE, scale="row", density.info="none",trace="none",
          key=TRUE, margins=c(10,10))
dev.off()

#gustatory receptors
#par(mar=c(4.1,4.1,4.1,4.1))
#png(file="heatmap_TCO_gust_receptors.png",height=2800,width=2800)
#gustatory_list <- as.character(rownames(gustatory_table))
#promoter_list <- as.character(rownames(top_table_e))
#gustatory_rows <- match(gustatory_list, promoter_list)
#length(gustatory_rows)
#head(gustatory_rows)
#gustatory_rows <- na.omit(gustatory_rows)
#selected  <- rownames(top_table_e[gustatory_rows])
#esetSel <- dp_eset[selected, ]
#heatmap.2(exprs(esetSel), symm=FALSE,symkey=FALSE, scale="row", density.info="none",trace="none",
#          key=TRUE, margins=c(10,10))
