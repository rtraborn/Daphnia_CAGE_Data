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
Dp_dge <- estimateTagwiseDisp(Dp_dge, trend="none")

plotBCV(Dp_dge)

et <- exactTest(Dp_dge)
topTags(et, n=20)
de_pE_male <- decideTestsDGE(et, adjust="BH",p.value=0.01)
summary(de_pE_male)
detags <- rownames(Dp_dge)[as.logical(de_pE_male)]
de.genes_male <- rownames(Dp_dge)[as.logical(de_pE_male)]
plotSmear(Dp_dge, de.tags = de.genes_male, cex = 0.5)
abline(h = c(-2, 2), col = "blue")

top_table_e <- topTags(et, p.value=0.01)

v <- voom(Dp_dge,design,plot=TRUE)
fit <- lmFit(v,design=design)
fit2 <- eBayes(fit)

topTable(fit2,coef=ncol(design))

#volcanoplot(fit2)

options(digits=3)
top_table <- topTable(fit2,coef=ncol(design),n=Inf,sort.by="p",adjust="BH",p=0.01)
results <- decideTests(fit2)
write.table(results,file="de_TCO_decideTest.txt",row.names=TRUE,quote=FALSE)

vennDiagram(results, names=c("Asexual females","Sexual females","Sexual Males"),include="up",circle.col=c("green","red","blue"))

#top_table2 <- cbind(top_table,results)
#write.table(top_table2,file="de_TCO_combined_table.txt", row.names=TRUE,quote=FALSE)

#number of differentially-related promoters
sum(top_table_e$adj.P.Val<0.01)

de_data <- Dp_dge$pseudo.counts

#differential analysis results
de_data <- cbind(de_data, et$table)

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

write.table(de_data,file="de_data_TCO_pE_male.txt",col.names=TRUE,quote=FALSE)

###########################################################################
#Adding gene annotation to promoters

#first, creating a genomicRanges object from the promoter data
de_GR <- with(de_data, GRanges(chr,
                                    IRanges(start, end, names=row.names(de_data)), strand))

#what does the GR object look like?
de_GR

#importing the gene annotation file
dpulex_genes <- read.table(file="/home/rtraborn/Daphnia/Daphnia_CAGE_Data/gene_annotations/dpulex_genes.bed", header=FALSE)
names(dpulex_genes) <- c("chr","start","end","geneID","score","strand","version","type","placeholder","ID2")

#region around gene to call a TSS
span <- 200

#adding space to the genes on the positive strand
dpulex_genes[dpulex_genes$strand=='+','start']  <- dpulex_genes[dpulex_genes$strand=='+',"start"]-span
dpulex_genes[dpulex_genes$strand=='+','end']  <- dpulex_genes[dpulex_genes$strand=='+',"end"]+span

#adding space to the genes on the negative strand
dpulex_genes[dpulex_genes$strand=='-','start']  <- dpulex_genes[dpulex_genes$strand=='-',"start"]-span
dpulex_genes[dpulex_genes$strand=='-','end']  <- dpulex_genes[dpulex_genes$strand=='-',"end"]+span

#create GRanges object for all genes
genes_GR <- with(dpulex_genes, GRanges(chr, IRanges(start, end, names=geneID),strand))

#overlap TSRs with all tag clusters
tc_overlaps <- findOverlaps(de_GR, genes_GR)

#store number of overlaps
tc_overlaps_count <- countOverlaps(de_GR, genes_GR)

match_hit2 <- as.data.frame(tc_overlaps)

#name the columns
names(match_hit2) <- c('query','subject')

promoter_index <- match_hit2$query
gene_index <- match_hit2$subject

#gene_names <- dpulex_genes[gene_index,'geneID']
promoter_table <- data.frame(de_data[promoter_index,])
promoter_IDs <- dpulex_genes[gene_index,"geneID"]
rownames(promoter_table) = make.names(promoter_IDs, unique=TRUE)

#remove duplicated entries
#match_hit2 <- match_hit2[!duplicated(match_hit2$query),]

write.table(promoter_table,file="TCO_promoter_de_table.txt",col.names=TRUE, row.names=TRUE)

###########################################################################
#Making heatmaps from the eset data we've generated

#par(mar=c(2.1,4.1,2.1,4.1))
#png(file="heatmap_TCO_all.png",height=1600,width=1200)
#selected  <- which(de_data$p.value[,2] <0.01)
#esetSel <- dp_eset[selected, ]
#heatmap(exprs(esetSel))
#dev.off()

par(mar=c(2.1,4.1,2.1,4.1))
png(file="heatmap_TCO_upreg1.png",height=1600,width=1200)
selected  <- which(de_pE_male[,1]==1)
esetSel <- dp_eset[selected, ]
heatmap(exprs(esetSel))
dev.off()
