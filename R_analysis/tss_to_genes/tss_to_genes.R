#Associates Tag Clusters (TSRs/promoter) to genes

#load libraries
library(GenomicRanges)
library(GenomicFeatures)

dpulex_genes <- read.table(file="/home/rtraborn/Daphnia/Daphnia_CAGE_Data/gene_annotations/dpulex_genes.bed", header=FALSE)
dpulex_promoters <- read.table(file="/home/rtraborn/Daphnia/CAGE/CAGEr_promoter_annotations/All.samples.consensusClusters.qLow0.1_qUp0.9.bed",header=FALSE,skip=1,fill=TRUE)

names(dpulex_genes) <- c("chr","start","end","geneID","score","strand","version","type","placeholder","ID2")

names(dpulex_promoters) <- c("chr","start","end","placeholder","score","strand","start2","end2","score2","score3","score4","score5")

#region around gene to call a TSS
span <- 500

print(head(dpulex_genes))

print(head(dpulex_promoters))


#adding space to the genes on the positive strand
dpulex_genes[dpulex_genes$strand=='+',"start"]  <- dpulex_genes[dpulex_genes$strand=='+',"start"]-span
dpulex_genes[dpulex_genes$strand=='+',"end"]  <- dpulex_genes[dpulex_genes$strand=='+',"end"]+span

#adding space to the genes on the negative strand
dpulex_genes[dpulex_genes$strand=='-',"start"]  <- dpulex_genes[dpulex_genes$strand=='-',"start"]-span
dpulex_genes[dpulex_genes$strand=='-',"end"]  <- dpulex_genes[dpulex_genes$strand=='-',"end"]+span

#create GRanges object for all tag clusters
promoters_GR <- with(dpulex_promoters, GRanges(chr, IRanges(as.numeric(start),as.numeric(end),names=placeholder,strand)))

genes_GR <- with(dpulex_genes, GRanges(chr, IRanges(start, end, names=geneID,strand)))

#overlap TSRs with all tag clusters
tc_overlaps <- findOverlaps(promoters_GR, genes_GR)

#store number of overlaps
tc_overlaps_count <- countOverlaps(promoters_GR, genes_GR)

#store number of overlaps in the big table of results
#data2$oc <- as.vector(tc_overlaps_count)

#create data frame that matches the Entrez TSS to the tag cluster
match_hit2 <- data.frame(names(promoters_GR)[queryHits(tc_overlaps)],
                        names(genes_GR)[subjectHits(tc_overlaps)],
                        stringsAsFactors=F
                        )

#name the columns
names(match_hit2) <- c('query','subject')

#remove duplicated entries
match_hit2 <- match_hit2[!duplicated(match_hit2$query),]

#name the rows
row.names(match_hit2) <- match_hit2$query

write.table(match_hit2,file="gene_names.txt",col.names=TRUE, row.names=TRUE)

q()

