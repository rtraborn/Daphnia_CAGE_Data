
#creates MA plot for asexual females and sexual females

require("CAGEr")

require("limma")

require("edgeR")

require("Biobase")

require("gplots")

require(CAGEr)

require(GenomicRanges)

require(GenomicFeatures)

dpulex_genes <- read.table(file="/home/rtraborn/Daphnia/Daphnia_CAGE_Data/gene_annotations/dpulex_genes.bed", header=FALSE)

names(dpulex_genes) <- c("chr","start","end","geneID","score","strand","version","type","placeholder","ID2")

head(dpulex_genes)

#importing the Dp CAGEr object

load("/home/rtraborn/Daphnia/Daphnia_CAGE_Data/R_analysis/promoter_calling_pipelines/TCO/Dp_CAGEset_TCO.RData")

#creating a count table

Dp_edger <- myCAGEset@consensusClustersTpmMatrix

#creating identifier for each tag cluster

Dp_edger_consensus_cluster <- consensusClusters(myCAGEset)

head(Dp_edger_consensus_cluster)

dim(Dp_edger_consensus_cluster)

rownames(Dp_edger) <- paste(Dp_edger_consensus_cluster$chr,
                                        Dp_edger_consensus_cluster$start,
                                        Dp_edger_consensus_cluster$end,
                                        Dp_edger_consensus_cluster$strand,
                                        sep="_")

#region around gene to call a TSS

span <- 500

#adding space to the genes on the positive strand  

dpulex_genes[dpulex_genes$strand=='+','start']  <- dpulex_genes[dpulex_genes$strand=='+',"start"]-span

dpulex_genes[dpulex_genes$strand=='+','end']  <- dpulex_genes[dpulex_genes$strand=='+',"end"]+span

#adding space to the genes on the negative strand 

dpulex_genes[dpulex_genes$strand=='-','start']  <- dpulex_genes[dpulex_genes$strand=='-',"start"]-span

dpulex_genes[dpulex_genes$strand=='-','end']  <- dpulex_genes[dpulex_genes$strand=='-',"end"]+span

#create GRanges object for all tag clusters

promoters_GR <- with(Dp_edger_consensus_cluster, GRanges(chr, IRanges(start,end,names=rownames(Dp_edger)),strand))

promoters_GR

genes_GR <- with(dpulex_genes, GRanges(chr, IRanges(start, end, names=geneID),strand))

#overlap TSRs with all tag clusters 

tc_overlaps <- findOverlaps(promoters_GR, genes_GR)

#store number of overlaps 

tc_overlaps_count <- countOverlaps(promoters_GR, genes_GR)

match_hit2 <- as.data.frame(tc_overlaps)

#name the columns 

names(match_hit2) <- c('query','subject')

promoter_index <- match_hit2$query

gene_index <- match_hit2$subject

#gene_names <- dpulex_genes[gene_index,'geneID']

promoter_table <- data.frame(Dp_edger_consensus_cluster[promoter_index,1:3])

edger_table <- data.frame(Dp_edger[promoter_index,])

head(edger_table)

promoter_table$strand <- Dp_edger_consensus_cluster$strand[promoter_index]

promoter_table$gene <- dpulex_genes[gene_index,4]

promoter_table$gene_start <- dpulex_genes[gene_index,'start']

promoter_table$gene_end <- dpulex_genes[gene_index,'end']

promoter_table$gene_strand <- dpulex_genes[gene_index,'strand']

promoter_table$expression <- Dp_edger_consensus_cluster$tpm[promoter_index]

head(Dp_edger_consensus_cluster)

head(promoter_table)

#head(dp_eset)

dp_eset <-new("ExpressionSet", exprs=as.matrix(edger_table))

dim(Dp_edger)

dim(promoter_table)

Dp_edger[1:20,]

rownames(edger_table) <- rownames(Dp_edger[promoter_index])

lib_sizes <- librarySizes(myCAGEset)

design <- model.matrix(~0+factor(c(1,1,1,3,3,2,2,2)))

colnames(design) <- c("mat_fem", "pE_fem", "mat_male")

design

contrasts.matrix <- makeContrasts(
    MvMf=mat_male-mat_fem,
    pEvM=pE_fem-mat_male,
    MfvpE=mat_fem-pE_fem,
    MvF=mat_male-(mat_fem+pE_fem)/2,
    AvS=mat_fem-(mat_fem+pE_fem)/2,
    levels=design)

contrasts.matrix

group <- c(rep("mat_fem",3),rep("mat_male",2),rep("pE_fem",3))

p_cutoff <- 0.01

rowsum_threshold <- 40

Dp_dge <- DGEList(counts=Dp_edger,group=group)

head(Dp_dge)

A <- rowSums(Dp_dge$counts) 

Dp_dge <- Dp_dge[A>rowsum_threshold,] 

head(Dp_dge$counts)

Dp_dge <- calcNormFactors(Dp_dge)

Dp_dge <- estimateCommonDisp(Dp_dge, verbose=T)

Dp_dge <- estimateTagwiseDisp(Dp_dge, trend="none")

plotBCV(Dp_dge)

v <- voom(Dp_dge, design, plot=TRUE)

fit <- lmFit(v,design)

fit2 <- contrasts.fit(fit,contrasts.matrix)

fit2 <- eBayes(fit2)

res <- decideTests(fit2,p.value=0.001,lfc=log2(2))

ind = which(apply(res,1,function(x) {length(which(x != 0))>0}) == T)

length(ind)

isDE <- as.logical(res)

isDE

DEnames <- rownames(fit2)[isDE]

head(res)

options(digits=3)

de_table2 <- topTable(fit2, coef=2, sort="none",adjust.method="BH",n=Inf)

plotMDS(v, labels=c("mf1","mf2","mf3","m1","m2","pE1","pE2","pE3"), main="MDS plot for all eight libraries")

fit$genes <- rownames(fit$coefficients)

head(fit$genes)

fit2[1438,]

#?volcanoplot

is(fit2)

#?plotMA

Dp_dge$samples$group

et <- exactTest(Dp_dge,pair=c("mat_fem","pE_fem"))

#head(et)

summary(de <- decideTestsDGE(et, p=0.01, adjust="BH"))

detags <- rownames(Dp_dge)[as.logical(de)]

png(file="mat_fem_v_pE_plotsmear.png",res=300)

plot.new()

plotSmear(et, de.tags=detags, smooth.scatter = FALSE, smearWidth = 0.5, lowess=FALSE); abline(h = c(-2, 2), col = "dodgerblue", lwd = 2,lty=3)

dev.off()

