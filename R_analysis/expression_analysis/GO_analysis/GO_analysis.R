#does hypergeometic analysis using upregulated genes
#modeled on second half of analysis.R from Dave Tang's 'ccl' repo

library(topGO)

### importing the gene annotation master table
Dp_JGI_genes <- read.table(file="/home/rtraborn/Daphnia/Daphnia_CAGE_Data/gene_annotations/Dp_JGI_genes_GO.txt",header=TRUE,stringsAsFactors=FALSE,sep="\t",skip=1)

Dp_JGI_genes <- c("promoter_ID","chr","start","end","strand","gene","GO_ID","KOG_JGI","KOG_ID","FlyBase_id_representative_homolog","interprot_representative_homolog")

## importing the de data from limma

## Males vs Asexual females
MvMF_de <- read.table(file="/home/rtraborn/Daphnia/Daphnia_CAGE_Data/R_analysis/expression_analysis/de_tables/TCO_MvMf_de_table_genes.txt",header=TRUE,stringsAsFactors=FALSE)

cutoff <- 0.001

MvMF_de$FDR <- as.numeric(MvMF_de$FDR)
MvMF_de$logFC <- as.numeric(MvMF_de$FDR)

#tag clusters significantly upregulated in Males
MvMF_up_list <- subset(MvMF_de, FDR<cutoff & logFC>0)
up_MvMF_gene <- MvMF_up_list$gene

length(up_MvMF_gene)
head(MvMF_up_list)

#universal list
all_promoter <- MvMF_up_list$promoter_ID
geneID2GO <- readMappings(file = "/home/rtraborn/Daphnia/Daphnia_CAGE_Data/R_analysis/expression_analysis/GO_analysis/dp_jgi_v11.map")

geneNames <- names(geneID2GO)
geneList <- factor(as.integer(geneNames %in% up_MvMF_gene))
names(geneList) <- geneNames
str(geneList)

#making an S4 object from topGO

dpGOdata <- new("topGOdata",
description = "Genes strongly upregulated in males", ontology = "BP",
allGenes = geneList,
nodeSize = 10,
annot = annFUN.gene2GO,
gene2GO = geneID2GO)

resultFisher <- runTest(dpGOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(dpGOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(dpGOdata, algorithm = "elim", statistic = "ks")

allRes <- GenTable(dpGOdata, classicFisher = resultFisher,
classicKS = resultKS, elimKS = resultKS.elim,
 orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)

printGraph(dpGOdata, resultKS, firstSigNodes = 10, useInfo = "all", pdfSW = TRUE,fn.prefix="MvMF_upreg_tGO_BP")
write.table(allRes,file="MvMF_upreg_GO_fisher_BP.txt",col.names=TRUE,row.names=FALSE,quote=TRUE,sep="\t")

dpGOdata <- new("topGOdata",
description = "Genes strongly upregulated in males", ontology = "MF",
allGenes = geneList,
nodeSize = 10,
annot = annFUN.gene2GO,
gene2GO = geneID2GO)

resultFisher <- runTest(dpGOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(dpGOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(dpGOdata, algorithm = "elim", statistic = "ks")

allRes <- GenTable(dpGOdata, classicFisher = resultFisher,
classicKS = resultKS, elimKS = resultKS.elim,
 orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)

printGraph(dpGOdata, resultKS, firstSigNodes = 10, useInfo = "all", pdfSW = TRUE,fn.prefix="MvMF_upreg_tGO_MF")
write.table(allRes,file="MvMF_upreg_GO_fisher_MF.txt",col.names=TRUE,row.names=FALSE,quote=TRUE,sep="\t")


## pE vs Asexual females
pEvMF_de <- read.table(file="/home/rtraborn/Daphnia/Daphnia_CAGE_Data/R_analysis/expression_analysis/de_tables/TCO_pEvM_de_table_genes.txt",header=TRUE,stringsAsFactors=FALSE)

cutoff <- 0.001

pEvMF_de$FDR <- as.numeric(pEvMF_de$FDR)
pEvMF_de$logFC <- as.numeric(pEvMF_de$FDR)

#tag clusters significantly upregulated in pE females relative to asexuals
pEvMF_up_list <- subset(pEvMF_de, FDR<cutoff & logFC>0)
up_pEvMF_gene <- pEvMF_up_list$gene

length(up_pEvMF_gene)
head(pEvMF_up_list)

#universal list
all_promoter <- pEvMF_up_list$promoter_ID
geneID2GO <- readMappings(file = "/home/rtraborn/Daphnia/Daphnia_CAGE_Data/R_analysis/expression_analysis/GO_analysis/dp_jgi_v11.map")

geneNames <- names(geneID2GO)
geneList <- factor(as.integer(geneNames %in% up_pEvMF_gene))
names(geneList) <- geneNames
str(geneList)

#making an S4 object from topGO

dpGOdata <- new("topGOdata",
description = "Genes strongly upregulated in pE vs asexual females", ontology = "BP",
allGenes = geneList,
nodeSize = 10,
annot = annFUN.gene2GO,
gene2GO = geneID2GO)

resultFisher <- runTest(dpGOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(dpGOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(dpGOdata, algorithm = "elim", statistic = "ks")

allRes <- GenTable(dpGOdata, classicFisher = resultFisher,
classicKS = resultKS, elimKS = resultKS.elim,
 orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)

printGraph(dpGOdata, resultKS, firstSigNodes = 10, useInfo = "all", pdfSW = TRUE,fn.prefix="pEvMF_upreg_tGO_BP")
write.table(allRes,file="upreg_pE_v_matfem_GO_fisher_BP.txt",col.names=TRUE,row.names=FALSE,quote=TRUE,sep="\t")

dpGOdata <- new("topGOdata",
description = "Genes strongly upregulated in pE vs asexual females", ontology = "MF",
allGenes = geneList,
nodeSize = 10,
annot = annFUN.gene2GO,
gene2GO = geneID2GO)

resultFisher <- runTest(dpGOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(dpGOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(dpGOdata, algorithm = "elim", statistic = "ks")

allRes <- GenTable(dpGOdata, classicFisher = resultFisher,
classicKS = resultKS, elimKS = resultKS.elim,
 orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)

printGraph(dpGOdata, resultKS, firstSigNodes = 10, useInfo = "all", pdfSW = TRUE,fn.prefix="pEvMF_upreg_tGO_MF")
write.table(allRes,file="upreg_pE_v_matfem_GO_fisher_MF.txt",col.names=TRUE,row.names=FALSE,quote=TRUE,sep="\t")

############ upregulated in asexual vs sexual females

pEvMF_de <- read.table(file="/home/rtraborn/Daphnia/Daphnia_CAGE_Data/R_analysis/expression_analysis/de_tables/TCO_pEvM_de_table_genes.txt",header=TRUE,stringsAsFactors=FALSE)

#tag clusters significantly upregulated in asexuals 
#(down in pE means up in asexuals)
pEvMF_down_list <- subset(pEvMF_de, FDR<cutoff & logFC<0)
down_pEvMF_gene <- pEvMF_down_list$gene

length(down_pEvMF_gene)
head(pEvMF_down_list)

geneID2GO <- readMappings(file = "/home/rtraborn/Daphnia/Daphnia_CAGE_Data/R_analysis/expression_analysis/GO_analysis/dp_jgi_v11.map")

geneNames <- names(geneID2GO)
geneList <- factor(as.integer(geneNames %in% down_pEvMF_gene))
names(geneList) <- geneNames
str(geneList)

#making an S4 object from topGO

dpGOdata <- new("topGOdata",
description = "Genes strongly upregulated in asexual females vs sexual (pE) females", ontology = "BP",
allGenes = geneList,
nodeSize = 10,
annot = annFUN.gene2GO,
gene2GO = geneID2GO)

resultFisher <- runTest(dpGOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(dpGOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(dpGOdata, algorithm = "elim", statistic = "ks")

allRes <- GenTable(dpGOdata, classicFisher = resultFisher,
classicKS = resultKS, elimKS = resultKS.elim,
 orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)

printGraph(dpGOdata, resultKS, firstSigNodes = 10, useInfo = "all", pdfSW = TRUE,fn.prefix="MfvpE_upreg_tGO_BP")
write.table(allRes,file="upreg_matfem_v_pE_GO_fisher_BP.txt",col.names=TRUE,row.names=FALSE,quote=TRUE,sep="\t")

dpGOdata <- new("topGOdata",
description = "Genes strongly upregulated in asexual females vs sexual (pE) females", ontology = "MF",
allGenes = geneList,
nodeSize = 10,
annot = annFUN.gene2GO,
gene2GO = geneID2GO)

resultFisher <- runTest(dpGOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(dpGOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(dpGOdata, algorithm = "elim", statistic = "ks")

allRes <- GenTable(dpGOdata, classicFisher = resultFisher,
classicKS = resultKS, elimKS = resultKS.elim,
 orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)

printGraph(dpGOdata, resultKS, firstSigNodes = 10, useInfo = "all", pdfSW = TRUE,fn.prefix="MfvpE_upreg_tGO_MF")
write.table(allRes,file="upreg_matfem_v_pE_GO_fisher_MF.txt",col.names=TRUE,row.names=FALSE,quote=TRUE,sep="\t")
