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

Dp_dge <- DGEList(counts=Dp_edger)

Dp_dge <- calcNormFactors(Dp_dge)

v <- voom(Dp_dge,design,plot=TRUE)

fit <- lmFit(v,design)
fit <- eBayes(fit)
volcanoplot(fit)

options(digits=3)
top <- topTable(fit,coef=2,n=Inf,sort.by="P",adjust="BH")

#number of differentially-related promoters
sum(top$adj.P.Val<0.01)

write.table(top,file="Dp_top_DGE_list.txt",col.names=TRUE,row.names=TRUE)

#plotMDS(v, labels=1:8, col=as.numeric(genotype))
