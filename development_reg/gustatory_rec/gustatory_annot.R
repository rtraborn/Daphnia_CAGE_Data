require(GenomicRanges)

#importing the gene annotation file
dpulex_genes <- read.table(file="/home/rtraborn/Daphnia/Daphnia_CAGE_Data/gene_annotations/dpulex_genes.bed", header=FALSE)
names(dpulex_genes) <- c("chr","start","end","geneID","score","strand","version","type","placeholder","ID2")

#create GRanges object for all genes
genes_GR <- with(dpulex_genes, GRanges(chr, IRanges(start, end, names=geneID),strand))

#importing the gustatory receptor gene annotation file
gr_genes <- read.table(file="/home/rtraborn/Daphnia/Daphnia_CAGE_Data/development_reg/gustatory_rec/Dp_gr_table.bed", header=FALSE, skip=1)
names(gr_genes) <- c("chr","start","end","strand","gene_ID","protein_ID","symbol")

gr_genes$start <- as.numeric(gr_genes$start)
gr_genes$end <- as.numeric(gr_genes$end)

#create GRanges object for all genes
gust_GR <- with(gr_genes, GRanges(chr, IRanges(start, end, names=protein_ID)))

#overlap TSRs with all tag clusters
gust_overlaps <- findOverlaps(gust_GR, genes_GR)

match_hit2 <- as.data.frame(gust_overlaps)
names(match_hit2) <- c('query','subject')

#remove duplicated entries
#match_hit2 <- match_hit2[!duplicated(match_hit2$query),]
match_hit2 <- match_hit2[!duplicated(match_hit2$subject),]

gust_index <- match_hit2$query
gene_index <- match_hit2$subject

head(match_hit2)

#gene_names <- dpulex_genes[gene_index,'geneID']
gust_table <- data.frame(dpulex_genes[gene_index,])
rownames(gust_table) <- gust_table$symbol
write.table(gust_table, file="Dpulex_gustatory_receptor_genes.bed",col.names=TRUE,row.names=FALSE,quote=FALSE)
