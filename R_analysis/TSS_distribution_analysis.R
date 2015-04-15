#this workflow generates TSS distribution figures around genomic intervals

#installing packages; uncomment if needed
#install.packages( c("data.table","plyr","reshape2","ggplot2","gridBase","devtools"))
#source("http://bioconductor.org/biocLite.R")
#biocLite(c("GenomicRanges","rtracklayer","impute","Rsamtools"))

library("devtools")
library("impute")
library("GenomicRanges")

install.packages("/home/rtraborn/genome_analysis/genomation", repos = NULL, type="source")

library("genomation")

setwd("/DATA/GROUP/rtraborn/Daphnia_Project/Dp_data")

bam.files <- list.files(path="/DATA/GROUP/rtraborn/Daphnia_Project/Dp_data/merged_bamfiles/",pattern="bam$")
print(bam.files)

#importing Dpulex genes
Dp_genes <- "FrozenGeneCatalog.bed"
Dp_sar <- "Daphnia_pulex_Sar_annot.bed"
Dpulex_gr <- readTranscriptFeatures(Dp_genes,remove.unusual=FALSE,unique.prom=TRUE)
#names(Dpulex_flank_500) <- c("chr","start","end","id","blank","strand","blank2","feature","blank3","id2")
Dp_Sar_gr <- readTranscriptFeatures(Dp_sar,remove.unusual=FALSE,unique.prom=TRUE)

Dpulex_promoters <- read.table(file="pre-ephippial_females.tagClusters.qLow0.1_qUp0.9.bed",header=F,stringsAsFactors=F)
names(Dpulex_promoters) <- c("chr","start","end","blank","score","strand","start2","end2","score2","score3","ol","ol2")

#creating a GRanges object from the genes
Dpulex_promoters_gr <- with(Dpulex_promoters, GRanges(chr, IRanges(start, end), strand=strand))

Dpulex_promoters_gr

setwd("/DATA/GROUP/rtraborn/Daphnia_Project/Dp_data/merged_bamfiles")

#analyzing the intersection between the CAGE and annotation files
CAGE_annot <- annotateWithGeneParts(Dpulex_promoters_gr, Dpulex_gr, intersect.chr = TRUE)

CAGE_annot

Dp_Sar_annot <- annotateWithGeneParts(Dpulex_promoters_gr, Dp_Sar_gr, intersect.chr = TRUE)

Dp_Sar_annot

#creating the Score Matrix Object
#sm <- ScoreMatrix(target=bam.files[1], windows=Dp_flank_500_gr, type = "bam")
#sm.scaled <- scaleScoreMatrix(sm)

#creating the pileup heatmap figure
#png(file="Dp_TSS_flank_500_heatmap2.png",height=960,width=960)

#creating the heatmap
#heatMatrix(sm.scaled, xcoords = c(-500, 0),kmeans=T,k=3, col = c("lightgray", "blue"))

#dev.off()
q()
