#creates a series of plots relating to differential expression analysis
# uses the bioconductor package 'limma'

require("CAGEr")
require("limma")
require("edgeR")

#importing the Dp CAGEr object
#load("DpCAGE_rep.RData")

#importing the edgeR-adapted consensus cluster counts file
Dp_edger <- read.table(file="Dp_rep_norm_cons_clusters.txt", header=TRUE, stringsAsFactors=FALSE)

#importing the design object
load("Dp_design.RData")

#checking to see what the data.frame looks like
head(Dp_edger)

#lib_sizes <- librarySizes(DpCAGE_rep)
#contrast.matrix <- makeContrasts(mat_male-mat_fem, pE_fem-mat_male, pE_fem-mat_fem, levels=design)

group <- c(rep("mat_fem",3),rep("mat_male",2),rep("pE_fem",3)) 
p_cutoff <- 0.01
rowsum_threshold <- 20

Dp_dge <- DGEList(counts=Dp_edger,group=group)
A <- rowSums(Dp_dge$counts) 
Dp_dge <- Dp_dge[A>rowsum_threshold,]  
Dp_dge <- calcNormFactors(Dp_dge)
Dp_dge <- estimateCommonDisp(Dp_dge, verbose=T)

de.tgw2_pE_male <- exactTest(Dp_dge,pair=c("pE_fem","mat_male"))

v <- voom(Dp_dge,design,plot=TRUE)
fit <- lmFit(v,design)
fit <- eBayes(fit)

volcanoplot(fit)

options(digits=3)
top_pE <- topTable(fit,coef="pE_fem",n=Inf,sort.by="p",adjust="BH",p=0.01)
top_mat_fem <- topTable(fit,coef="mat_fem",n=Inf,sort.by="p",adjust="BH",p=0.01)
top_mat_male <- topTable(fit,coef="mat_male",n=Inf,sort.by="p",adjust="BH",p=0.01)

#number of differentially-related promoters
sum(top_pE$adj.P.Val<0.01)
sum(top_mat_fem$adj.P.Val<0.01)
sum(top_mat_male$adj.P.Val<0.01)

write.table(top_pE,file="Dp_top_pE.txt",col.names=TRUE,row.names=TRUE)
write.table(top_mat_fem,file="Dp_top_mat_fem.txt",col.names=TRUE,row.names=TRUE)
write.table(top_mat_male,file="Dp_de_male.txt",col.names=TRUE,row.names=TRUE)

de_data <- Dp_dge$pseudo.counts

#differential analysis results
de_data <- cbind(de_data, de.tgw2_pE_male$table)

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

#creating a genomicRanges object
de_data_gr <- with(de_data, GRanges(chr,
                                    IRanges(start, end, names=row.names(de_data)), strand))

#what does the GR object look like?
de_data_gr

q()
#annotation (in progress- needs to be converted to Daphnia; this is a human example)
#columns <- c("GeneID","Symbol", "chromosome", "type_of_gene")
#GeneInfo <- read.columns("130811-Homo_sapiens.gene_info",required=columns, stringsAsFactors = FALSE)
#m <- match(gene$annotation$GeneID, GeneInfo$GeneID)

#Ann <- cbind(gene$annotation[, c("GeneID", "Chr", "Length")], GeneInfo[m, c("Symbol", "type_of_gene")])
#Ann$Chr <- unlist(lapply(strsplit(Ann$Chr, ";"), function(x) paste(unique(x), collapse = "|")))
#Ann$Chr <- gsub("chr", "", Ann$Chr)
#head(Ann)

