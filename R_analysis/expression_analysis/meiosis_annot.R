require(GenomicRanges)

#importing the gene annotation file
dpulex_genes <- read.table(file="/home/rtraborn/Daphnia/Daphnia_CAGE_Data/gene_annotations/dpulex_genes.bed", header=FALSE)
names(dpulex_genes) <- c("chr","start","end","geneID","score","strand","version","type","placeholder","ID2")

#create GRanges object for all genes
genes_GR <- with(dpulex_genes, GRanges(chr, IRanges(start, end, names=geneID),strand))

#importing the meiosis gene annotation file
meiosis_genes <- read.table(file="/home/rtraborn/Daphnia/Daphnia_CAGE_Data/development_reg/meiosis/meiosis_genes_schurko.bed", header=TRUE)
names(meiosis_genes) <- c("chr","start","end","name","score","strand")

#create GRanges object for all genes
meiosis_GR <- with(meiosis_genes, GRanges(chr, IRanges(start, end, names=name)))

#overlap TSRs with all tag clusters
meiosis_overlaps <- findOverlaps(meiosis_GR, genes_GR)

match_hit2 <- as.data.frame(meiosis_overlaps)
names(match_hit2) <- c('query','subject')

#remove duplicated entries
#match_hit2 <- match_hit2[!duplicated(match_hit2$query),]
match_hit2 <- match_hit2[!duplicated(match_hit2$subject),]

meiosis_index <- match_hit2$query
gene_index <- match_hit2$subject

head(match_hit2)

#gene_names <- dpulex_genes[gene_index,'geneID']
meiosis_table <- data.frame(dpulex_genes[gene_index,])
rownames(meiosis_table) <- meiosis_table$geneID
#rownames(meiosis_table) <- proteinID
#rownames(meiosis_table) <- gene_IDs
write.table(meiosis_table, file="Dpulex_meiosis_genes.bed",col.names=TRUE,row.names=TRUE,quote=FALSE)
