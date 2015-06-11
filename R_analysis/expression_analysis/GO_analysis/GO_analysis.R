#does hypergeometic analysis using upregulated genes
#modeled on second half of analysis.R from Dave Tang's 'ccl' repo

library(GO.db)
library(GOstats)

pE_vs_male_de <- read.table(file="/home/rtraborn/Daphnia/Daphnia_CAGE_Data/R_analysis/expression_analysis/de_tables/TCO_promoter_de_pE_male_master_table.txt",header=FALSE,stringsAsFactors=FALSE)

cutoff <- 0.01

#tag clusters differentially expressed with the CCL2 condition
up_pE_list <- subset(pE_vs_male_de, FDR<cutoff & logFC>0))

up_pE_promoter <- up_pE_list$promoter_ID
up_pE_gene <- up_pE_list$gene
up_pE_gene <- up_pE_list$GO_ID

#universal list
all_promoter <- up_pE_list$promoter_ID
all_gene <- up_pE_list$gene
all_GO <- up_pE_list$GO_ID


## start here once I figure out my next steps re: making an orgDb package (or not)
#tag clusters to Entrez gene 

up_ccl2_entrez2 <- unique(up_ccl2_entrez2)
universe_entrez2 <- match_hit2[universe2,2]
universe_entrez2 <- unique(universe_entrez2)

#load library
library(org.Hs.eg.db)

#parameters for hypergeometric test
ccl2_bp_param <- new('GOHyperGParams',
              geneIds=up_ccl2_entrez2,
              universeGeneIds=universe_entrez2,
              ontology='BP',
              pvalueCutoff=cutoff,
              conditional=F,
              testDirection='over',
              annotation="org.Hs.eg.db"
             )

#perform the test
ccl2_bp <- hyperGTest(ccl2_bp_param)

#store results
ccl2_bp_result <- summary(ccl2_bp)

#FDR
ccl2_bp_result$adj_p_value <- p.adjust(ccl2_bp_result$Pvalue, method="BH", n=1910)

#repeat analysis but for bFGF libraries
up_bfgf <- row.names(subset(data2, FDR<cutoff & logFC<0 & oc>0))
up_bfgf_entrez <- match_hit2[up_bfgf,2]
up_bfgf_entrez <- unique(up_bfgf_entrez)

bfgf_bp_param <- new('GOHyperGParams',
              geneIds=up_bfgf_entrez,
              universeGeneIds=universe_entrez2,
              ontology='BP',
              pvalueCutoff=cutoff,
              conditional=F,
              testDirection='over',
              annotation="org.Hs.eg.db"
             )
bfgf_bp <- hyperGTest(bfgf_bp_param)
bfgf_bp_result <- summary(bfgf_bp)
bfgf_bp_result$adj_p_value <- p.adjust(bfgf_bp_result$Pvalue, method="BH", n=1191)

#For hypoxia related genes
#load libraries
library(org.Hs.eg.db)
library(GO.db)

#create gene ontology object
go_object <- as.list(org.Hs.egGO2EG)

#all GO terms related to hypoxia via AmiGO
#http://amigo.geneontology.org/cgi-bin/amigo/search.cgi?search_query=hypoxia&search_constraint=term&action=new-search
#find all Entrez genes associated with these GO IDs
a <- go_object['GO:0097411']
b <- go_object['GO:0001666']
c <- go_object['GO:0070483']
d <- go_object['GO:0071456']
e <- go_object['GO:1900037']
f <- go_object['GO:1900038']
g <- go_object['GO:1900039']
h <- go_object['GO:1990144']
i <- go_object['GO:1902071']
j <- go_object['GO:1902072']
k <- go_object['GO:1902073']
l <- go_object['GO:0061418']
m <- go_object['GO:0061428']
n <- go_object['GO:0061419']
o <- go_object['GO:2000777']

#store all
all <- list(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o)

#keep unique list of Entrez genes
all <- unique(unlist(all, use.names=F))

#store as data frame
all <- data.frame(gene=all)

#store only annoated tag clusters
data_annotated <- data2[!is.na(data2$gene),]

#results of hypoxia related genes
all_merge <- merge(x=all, y=data_annotated, by='gene')
