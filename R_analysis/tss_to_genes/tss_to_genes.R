#Associates Tag Clusters (TSRs/promoter) to genes

#load libraries
library(GenomicRanges)
library(GenomicFeatures)

dpulex_genes <- read.table(file="/home/rtraborn/Daphnia/Daphnia_CAGE_Data/gene_annotations/dpulex_genes.bed", header=FALSE)
dpulex_promoters <- read.table(file="/home/rtraborn/Daphnia/CAGE/CAGEr_promoter_annotations/All.samples.consensusClusters.qLow0.1_qUp0.9.bed",header=FALSE,sep="\t",stringsAsFactors=FALSE)

names(dpulex_genes) <- c("chr","start","end","geneID","score","strand","version","type","placeholder","ID2")

names(dpulex_promoters) <- c("chr","start","end","placeholder","score","strand","start2","end2","score2","score3","score4","score5")

#region around gene to call a TSS
span <- 500

#adding space to the genes on the positive strand
dpulex_genes[dpulex_genes$strand=='+','start']  <- dpulex_genes[dpulex_genes$strand=='+',"start"]-span
dpulex_genes[dpulex_genes$strand=='+','end']  <- dpulex_genes[dpulex_genes$strand=='+',"end"]+span

#adding space to the genes on the negative strand
dpulex_genes[dpulex_genes$strand=='-','start']  <- dpulex_genes[dpulex_genes$strand=='-',"start"]-span
dpulex_genes[dpulex_genes$strand=='-','end']  <- dpulex_genes[dpulex_genes$strand=='-',"end"]+span

#create GRanges object for all tag clusters
promoters_GR <- with(dpulex_promoters, GRanges(chr, IRanges(start,end,names=placeholder),strand))

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
promoter_table <- data.frame(dpulex_promoters[promoter_index,1:3])
promoter_table$strand <- dpulex_promoters[promoter_index,6]

promoter_table$gene <- dpulex_genes[gene_index,4]

promoter_table$gene_start <- dpulex_genes[gene_index,'start']
promoter_table$gene_end <- dpulex_genes[gene_index,'end']
promoter_table$gene_strand <- dpulex_genes[gene_index,'strand']

#names(combined.table2) <- c("chr","start","end","gene_ID","strand","start_TSR","end_TSR")

#remove duplicated entries
#match_hit2 <- match_hit2[!duplicated(match_hit2$query),]

write.table(promoter_table,file="TCO_promoter_table.txt",col.names=TRUE, row.names=TRUE)

q()

