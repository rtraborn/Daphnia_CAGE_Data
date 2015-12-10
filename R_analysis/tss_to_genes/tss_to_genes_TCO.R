#Associates Tag Clusters (TSRs/promoter) to genes
#Modified to perform the conversion from the pooled tagClusters file: /home/rtraborn/Daphnia/Daphnia_CAGE_Data/R_analysis/promoter_calling_pipelines/TCO/tagClusters/pooled_samples
#Included "no gene"


#load libraries
library(GenomicRanges)
library(GenomicFeatures)

dpulex_genes <- read.table(file="/home/rtraborn/Daphnia/Daphnia_CAGE_Data/gene_annotations/dpulex_genes.bed", header=FALSE)
dpulex_promoters <- read.table(file="/home/rtraborn/Daphnia/Daphnia_CAGE_Data/R_analysis/promoter_calling_pipelines/TCO/tagClusters/pooled_samples/TSR_breadth_SI_9_2_15.txt", header=TRUE)

#dpulex_promoters <- read.table(file="/home/rtraborn/Daphnia/Daphnia_CAGE_Data/promoter_arch/D_pulex/TCO/ConsensusClusters_TCO.bed",header=FALSE,sep="\t",stringsAsFactors=FALSE,skip=1)
#dpulex_expression <- read.table(file="/home/rtraborn/Daphnia/Daphnia_CAGE_Data/promoter_arch/D_pulex/TCO/consensus_clusters_TCO.txt",stringsAsFactors=FALSE,header=TRUE) #uncomment this if you want to include this file (the pooled tagClusters file already has it)

names(dpulex_genes) <- c("chr","start","end","geneID","score","strand","version","type","placeholder","ID2")
names(dpulex_promoters) <- c("chr","start","end","breadth","placeholder","strand","nTSSs","SI")

#region around gene to call a TSS
span <- 500

#adding space to the genes on the positive strand
dpulex_genes[dpulex_genes$strand=='+','start']  <- dpulex_genes[dpulex_genes$strand=='+',"start"]-span
dpulex_genes[dpulex_genes$strand=='+','end']  <- dpulex_genes[dpulex_genes$strand=='+',"end"]

#adding space to the genes on the negative strand
dpulex_genes[dpulex_genes$strand=='-','start']  <- dpulex_genes[dpulex_genes$strand=='-',"start"]
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
promoter_table$breadth <- dpulex_promoters[promoter_index,"breadth"]
promoter_table$gene_strand <- dpulex_genes[gene_index,'strand']
promoter_table$expression <- dpulex_promoters[promoter_index,7]
promoter_table$SI <- dpulex_promoters[promoter_index,8]

promoter_table$promoter <- paste(dpulex_promoters[promoter_index,1], dpulex_promoters[promoter_index,2], dpulex_promoters[promoter_index,3], dpulex_promoters[promoter_index,6],sep="_")

rownames(promoter_table) <- 1:nrow(promoter_table)

#table for those promoters that don't match with genes
promoter_table2 <- data.frame(dpulex_promoters[-promoter_index,1:3])
promoter_table2$strand <- dpulex_promoters[-promoter_index,6]
promoter_table2$gene <- "No_gene"

promoter_table2$gene_start <- NA
promoter_table2$gene_end <- NA
promoter_table2$gene_strand <- NA
promoter_table2$breadth <- dpulex_promoters[-promoter_index,"breadth"]
promoter_table2$expression <- dpulex_promoters[-promoter_index,7]
promoter_table2$SI <- dpulex_promoters[-promoter_index,8]

promoter_table2$promoter <- paste(dpulex_promoters[-promoter_index,1], dpulex_promoters[-promoter_index,2], dpulex_promoters[-promoter_index,3], dpulex_promoters[-promoter_index,6],sep="_")

rownames(promoter_table2) <- 1:nrow(promoter_table2)

#names(combined.table2) <- c("chr","start","end","gene_ID","strand","start_TSR","end_TSR")
#remove duplicated entries
#match_hit2 <- match_hit2[!duplicated(match_hit2$query),]

#combining both tables
promoter_table_combined <- rbind(promoter_table, promoter_table2)
rownames(promoter_table_combined) <- 1:nrow(promoter_table_combined)

write.table(promoter_table,file="TC_list_table_genes.txt",col.names=TRUE, row.names=TRUE,quote=FALSE)
write.table(promoter_table2, file="TC_list_table_no_genes.txt", col.names=TRUE, row.names=TRUE, quote=FALSE)
write.table(promoter_table_combined, file="TC_list_table_complete.txt", col.names=TRUE, row.names=TRUE, quote=FALSE)

q()

