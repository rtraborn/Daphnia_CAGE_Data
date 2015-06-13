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

head(dp_eset)

#checking to see what the data.frame looks like
head(Dp_edger)

lib_sizes <- librarySizes(myCAGEset)
#contrast.matrix <- makeContrasts(mat_male-mat_fem, pE_fem-mat_male, pE_fem-mat_fem, levels=design)

design <- model.matrix(~ 0+factor(c(1,1,1,3,3,2,2,2)))
colnames(design) <- c("mat_fem", "pE_fem", "mat_male")
design

group <- c(rep("mat_fem",3),rep("mat_male",2),rep("pE_fem",3))
group
p_cutoff <- 0.01
rowsum_threshold <- 20

Dp_dge <- DGEList(counts=Dp_edger,group=group)
head(Dp_dge)
A <- rowSums(Dp_dge$counts) 
Dp_dge <- Dp_dge[A>rowsum_threshold,]  
#Dp_dge <- calcNormFactors(Dp_dge)
Dp_dge <- estimateCommonDisp(Dp_dge, verbose=T)
Dp_dge <- estimateTagwiseDisp(Dp_dge, trend="none")
plotBCV(Dp_dge)
v <- voom(Dp_dge, design, plot=TRUE)
fit <- lmFit(v,design)
fit <- eBayes(fit)
options(digits=3)
de_table <- topTable(fit, coef=ncol(design), sort.by="logFC",adjust.method="BH")
head(de_table)
plotMDS(v, labels=c("mfem1","mfem2","mfem3","male1","male2","pE1","pE2","pE3"), main="MDS plot for all eight libraries")
volcanoplot(fit,coef=ncol(design))

###### male vs pE females #######################
et <- exactTest(Dp_dge,pair=c("mat_male","pE_fem"))
topTags(et, n=20)
de_pE_male <- decideTestsDGE(et, adjust="BH",p.value=0.01)
summary(de_pE_male)
detags <- rownames(Dp_dge)[as.logical(de_pE_male)]
de.genes_pE_male <- rownames(Dp_dge)[as.logical(de_pE_male)]
plotSmear(Dp_dge, de.tags = de.genes_pE_male, cex = 0.5)
abline(h = c(-2, 2), col = "blue")

top_table_e <- topTags(et, sort.by="PValue",n=Inf)
top_table_e$FDR <- p.adjust(top_table_e$PValue, method = "BH")
dim(top_table_e)

data_coord2 <- matrix(data=unlist(strsplit(rownames(top_table_e), split="_")),
                      nrow= length(row.names(top_table_e)),
                      byrow=T)

data_coord2 <- as.data.frame(data_coord2, stringsAsFactors=F)
head(data_coord2)

col_1 <- data_coord2[,1]
col_2 <- data_coord2[,2]
chr_col <- paste(col_1,col_2,sep="_")
data_coord2 <- data_coord2[,-2]
data_coord2[,1] <- chr_col
data_coord2[,2] <- as.numeric(data_coord2[,2])
data_coord2[,3] <- as.numeric(data_coord2[,3])
names(data_coord2) <- c('chr','start','end','strand')

#first, creating a genomicRanges object from the promoter data
de_GR_1 <- with(data_coord2, GRanges(chr,
                                    IRanges(start, end, names=row.names(top_table_e)), strand))


de_GR_1 <- promoters(de_GR_1, upstream=200, downstream=200)

#what does the GR object look like?
de_GR_1

options(digits=3)

###### pE female vs mature females #######################
et2 <- exactTest(Dp_dge,pair=c("pE_fem","mat_fem"))
topTags(et2, n=20)
de_pE_matfem <- decideTestsDGE(et2, adjust="BH",p.value=0.01)
summary(de_pE_matfem)
detags <- rownames(Dp_dge)[as.logical(de_pE_matfem)]
de.genes_pE_matfem <- rownames(Dp_dge)[as.logical(de_pE_matfem)]
plotSmear(Dp_dge, de.tags = de.genes_pE_matfem, cex = 0.5)
abline(h = c(-2, 2), col = "blue")

top_table_e2 <- topTags(et2, sort.by="PValue",n=Inf)
top_table_e2$FDR <- p.adjust(top_table_e2$PValue, method = "BH")

data_coord2 <- matrix(data=unlist(strsplit(rownames(top_table_e2), split="_")),
                      nrow= length(row.names(top_table_e2)),
                      byrow=T)
data_coord2 <- as.data.frame(data_coord2, stringsAsFactors=F)
head(data_coord2)

col_1 <- data_coord2[,1]
col_2 <- data_coord2[,2]
chr_col <- paste(col_1,col_2,sep="_")
data_coord2 <- data_coord2[,-2]
data_coord2[,1] <- chr_col
data_coord2[,2] <- as.numeric(data_coord2[,2])
data_coord2[,3] <- as.numeric(data_coord2[,3])
names(data_coord2) <- c('chr','start','end','strand')

#first, creating a genomicRanges object from the promoter data
de_GR_2 <- with(data_coord2, GRanges(chr,
                                    IRanges(start, end, names=row.names(top_table_e2)), strand))

de_GR_2 <- promoters(de_GR_2, upstream=200, downstream=200)

#what does the GR object look like?
de_GR_2

###### asexual females vs males  #######################
et3 <- exactTest(Dp_dge,pair=c("mat_fem","mat_male"))
topTags(et3, n=20)
de_matfem_male <- decideTestsDGE(et3, adjust="BH",p.value=0.01)
summary(de_matfem_male)
detags <- rownames(Dp_dge)[as.logical(de_matfem_male)]
de.genes_matfem_male <- rownames(Dp_dge)[as.logical(de_matfem_male)]
plotSmear(Dp_dge, de.tags = de.genes_matfem_male, cex = 0.5)
abline(h = c(-2, 2), col = "blue")

top_table_e3 <- topTags(et3, sort.by="PValue",n=Inf)
top_table_e3$FDR <- p.adjust(top_table_e3$PValue, method = "BH")

data_coord2 <- matrix(data=unlist(strsplit(rownames(top_table_e3), split="_")),
                      nrow= length(row.names(top_table_e3)),
                      byrow=T)

data_coord2 <- as.data.frame(data_coord2, stringsAsFactors=F)
head(data_coord2)

col_1 <- data_coord2[,1]
col_2 <- data_coord2[,2]
chr_col <- paste(col_1,col_2,sep="_")
data_coord2 <- data_coord2[,-2]
data_coord2[,1] <- chr_col
data_coord2[,2] <- as.numeric(data_coord2[,2])
data_coord2[,3] <- as.numeric(data_coord2[,3])
names(data_coord2) <- c('chr','start','end','strand')

#first, creating a genomicRanges object from the promoter data
de_GR_3 <- with(data_coord2, GRanges(chr,
                                    IRanges(start, end, names=row.names(top_table_e3)), strand))


de_GR_3 <- promoters(de_GR_3, upstream=200, downstream=200)

#what does the GR object look like?
de_GR_3

#removing the limma commands for now
v <- voom(Dp_dge,design,plot=TRUE)
fit <- lmFit(v,design=design)
head(fit)
fit2 <- eBayes(fit)

topTable(fit2,coef=ncol(design))

top_table <- topTable(fit2,coef=ncol(design),n=Inf,sort.by="p",adjust="BH",p=0.01)
results <- decideTests(fit2)

vennDiagram(results, names=c("Asexual females","Sexual females","Sexual Males"),include="up",circle.col=c("green","red","blue"))

#top_table2 <- cbind(top_table,results)
#write.table(top_table2,file="de_TCO_combined_table.txt", row.names=TRUE,quote=FALSE)

head(et)

head(top_table_e)

################################### male vs pE females ############################

de_data <- Dp_dge$pseudo.counts

#differential analysis results
de_data <- cbind(de_data, et$table)

#calculating the false discovery rate
de_data$FDR <- p.adjust(de_data$PValue, method = 'BH')

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

#write.table(de_data1,file="TCO_pE_v_male_de.txt",col.names=TRUE,quote=FALSE)

################################### pE females vs asexual females ############################

de_data <- Dp_dge$pseudo.counts

#differential analysis results
de_data <- cbind(de_data, et2$table)

#calculating the FDR
de_data$FDR <- p.adjust(de_data$PValue, method = 'BH')

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

#write.table(de_data2,file="TCO_pE_v_matfem_de.txt",col.names=TRUE,quote=FALSE)

###################################  asexual females vs males  ############################

de_data <- Dp_dge$pseudo.counts

#differential analysis results
de_data <- cbind(de_data, et3$table)

#calculating the FDR
de_data$FDR <- p.adjust(de_data$PValue, method = 'BH')

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

#write.table(de_data3,file="TCO_matfem_v_males_de.txt",col.names=TRUE,quote=FALSE)

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
promoter_table <- data.frame(de_data1[promoter_index,])
gene_IDs <- dpulex_genes[gene_index,"geneID"]
promoter_table$gene <- gene_IDs

#remove duplicated entries
#match_hit2 <- match_hit2[!duplicated(match_hit2$query),]

#write.table(promoter_table,file="TCO_promoter_de_pE_male.txt",col.names=TRUE, row.names=TRUE)

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

#write.table(promoter_table,file="TCO_promoter_de_pE_matfem.txt",col.names=TRUE, row.names=TRUE)

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

#write.table(promoter_table,file="TCO_promoter_de_matfem_v_male.txt",col.names=TRUE, row.names=TRUE)

####################################

#number of differentially-related promoters
sum(de_data1$FDR<0.01)
sum(de_data2$FDR<0.01)
sum(de_data3$FDR<0.01)
sum(de_table$B<0.01)

################################# Gene Family Analyses ##############################

#importing the meiosis gene annotation file
meiosis_genes <- read.table(file="/home/rtraborn/Daphnia/Daphnia_CAGE_Data/development_reg/meiosis/Dpulex_meiosis_genes.bed", header=TRUE)
meiosis_IDs <- meiosis_genes$geneID
promoter_list <- promoter_table$gene
my_index <-  match(meiosis_IDs,promoter_table$gene)
length(my_index)
meiosis_table <- promoter_table[my_index,]
meiosis_table <- na.omit(meiosis_table)

#differentially-expressed meiosis genes only
de_index <- which(meiosis_table$de == 1)
meiosis_de <- meiosis_table[de_index,]

#head(meiosis_table)
#write.table(meiosis_table,file="meiosis_table.txt",col.names=TRUE,row.names=TRUE,quote=FALSE)
#write.table(meiosis_de,file="meiosis_table_de.txt",col.names=TRUE,row.names=TRUE,quote=FALSE)

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
selected  <- rownames(top_table_e)
esetSel <- dp_eset[selected,]
#heatmap.2(exprs(esetSel), symm=FALSE,symkey=FALSE,scale="row", density.info="none",trace="none",
#          key=TRUE,margins=c(10,10))
dev.off()

#male vs pE
par(mar=c(4.1,4.1,4.1,4.1))
png(file="heatmap_TCO_pE_v_male.png",height=2800,width=2800)
selected  <- rownames(top_table_e[1:200])
is(selected)
esetSel <- dp_eset[selected, ]
#heatmap.2(exprs(esetSel), symm=FALSE,symkey=FALSE,scale="row", density.info="none",trace="none",
#          key=TRUE,margins=c(10,10))
#dev.off()

#mat_fem vs pE
par(mar=c(4.1,4.1,4.1,4.1))
png(file="heatmap_TCO_matfem_v_pE.png",height=2800,width=2800)
selected  <- rownames(top_table_e2[1:200])
esetSel <- dp_eset[selected, ]
#heatmap.2(exprs(esetSel), symm=FALSE,symkey=FALSE,scale="row", density.info="none",trace="none",
#          key=TRUE,margins=c(10,10))
#dev.off()

#male vs mat_fem
par(mar=c(4.1,4.1,4.1,4.1))
png(file="heatmap_TCO_male_v_matfem.png",height=2800,width=2800)
selected  <- rownames(top_table_e3[1:200])
print(head(selected))
esetSel <- dp_eset[selected, ]
#heatmap.2(exprs(esetSel), symm=FALSE,symkey=FALSE, scale="row", density.info="none",trace="none",
#          key=TRUE, margins=c(10,10))
#dev.off()

#meiosis genes
par(mar=c(4.1,4.1,4.1,4.1))
png(file="heatmap_TCO_meiosis.png",height=2800,width=2800)
meiosis_rows <- match(rownames(meiosis_table), rownames(top_table_e))
length(meiosis_rows)
head(meiosis_rows)
meiosis_rows <- na.omit(meiosis_rows)
selected  <- rownames(top_table_e[meiosis_rows])
esetSel <- dp_eset[selected, ]
#heatmap.2(exprs(esetSel), symm=FALSE,symkey=FALSE, scale="row", density.info="none",trace="none",
#          key=TRUE, margins=c(10,10))
#dev.off()

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
