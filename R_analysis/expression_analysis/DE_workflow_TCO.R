# Analysis of Daphnia CAGE data, including promoter annotation using CAGEr and differential expression analysis
# Credit to Dave Tang for his github tutorial (located here: https://github.com/davetang/ccl2/blob/master/analysis.R)

#Session Info

sessionInfo()

require(CAGEr)

#loads D. pulex genome BS object
require(BSgenome.Dpulex.JGI.dpulex)

load("/home/rtraborn/Daphnia/Daphnia_CAGE_Data/R_analysis/promoter_calling_pipelines/Dp_TCO.RData")

getCTSS(myCAGEset, mappingQualityThreshold=10, sequencingQualityThreshold=-1)
cage_bam_ctss <- CTSStagCount(myCAGEset)

#tpm normalisation
normalizeTagCount(myCAGEset, method = "simpleTpm")

#tag clustering with 8 cores
clusterCTSS(object = myCAGEset,
            threshold = 1,
            thresholdIsTpm = TRUE,
            nrPassThreshold = 1,
            method = "distclu",
            maxDist = 20,
            removeSingletons = TRUE,
            keepSingletonsAbove = 5,
            useMulticore = T,
            nrCores = 8
           )

#cumulative sum of CAGE signal along genomic region
cumulativeCTSSdistribution(myCAGEset, clusters = "tagClusters", useMulticore = T, nrCores = 6)
quantilePositions(myCAGEset, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)

#aggregating tag clusters across multiple CAGE datasets
aggregateTagClusters(myCAGEset, tpmThreshold = 5, qLow = 0.1, qUp = 0.9, maxDist = 100)

#Extracting consensus clusters from CAGEset object
cage_bam_consensus_cluster <- consensusClusters(myCAGE_set)

#CAGE data based expression clustering
getExpressionProfiles(myCAGEset, what = "consensusClusters",
                      tpmThreshold = 10, nrPassThreshold = 1,
                      method = "som", xDim = 4, yDim = 2
                      )

#creating a new object for expression analysis
for_edger <- myCAGEset

for_edger_count <- for_edger@consensusClustersTpmMatrix

for_edger_consensus_cluster <- consensusClusters(for_edger) 

#promoter_table <- read.table(file="/home/rtraborn/Daphnia/Daphnia_CAGE_Data/data_tables/TCO/TCO_combined_cons_table.txt",header=TRUE,stringsAsFactors=FALSE)
           
for_edger_consensus_cluster$id <- paste(for_edger_consensus_cluster$chr,
					for_edger_consensus_cluster$start,
					for_edger_consensus_cluster$end,
					for_edger_consensus_cluster$strand,
                                        sep="_")

#set threshold for independent filtering
rowsum_threshold <- 20

#set p-value threshold for all analyses
cutoff <- 0.01

#preparing object for edgeR
d2 <- for_edger_count
rownames(d2) <- for_edger_consensus_cluster$id

#independent filtering
d2 <- d2[rowSums(d2)>rowsum_threshold,]

#load library
library(edgeR)

#Expression comparisons
group <- c(rep("asexual_females",3),rep("mature_males",2),rep("pE_females",3))
print(group)
print(ncol(d2))

d2 <- DGEList(counts = d2, group=group)

#TMM normalisation
d2 <- calcNormFactors(d2)

#calculate dispersion
d2 <- estimateCommonDisp(d2, verbose=T)
d2 <- estimateTagwiseDisp(d2)

#differential expression test
de.tgw2 <- exactTest(d2)

#summary of differentially expressed tag clusters
summary(decideTestsDGE(de.tgw2, p.value=cutoff))

#create data frame for results
#store normalised counts
data2 <- d2$pseudo.counts

#store differential analysis results
data2 <- cbind(data2, de.tgw2$table)

#calculate FDR
data2$FDR <- p.adjust(data2$PValue, method='BH')

#store dispersion of each tag cluster
data2$tw_dis <- d2$tagwise.dispersion

#prepare coordinates of each tag cluster
data_coord2 <- matrix(data=unlist(strsplit(row.names(data2), split="_")),
                      nrow= length(row.names(data2)),
                      byrow=T)
data_coord2 <- as.data.frame(data_coord2, stringsAsFactors=F)
names(data_coord2) <- c('chr','start','end','strand')

#store coordinates of each tag cluster
data2 <- cbind(data2, data_coord2)

#create column for differential expression status
#1 for DE and 0 for not
data2$de <- as.numeric(data2$FDR<cutoff)

#convert coordinates to numeric
data2$start <- as.numeric(data2$start)
data2$end <- as.numeric(data2$end)

write.table(data2,file="TCO_edgeR_dE_table.txt",quote=FALSE,row.names=TRUE,col.names=TRUE)
           
#load libraries
#library(GenomicRanges)
#library(GenomicFeatures)

#create GRanges object for all tag clusters
#data_grange2 <- with(data2,
#                    GRanges(chr, IRanges(start,
#                                         end,
#                                         names=row.names(data2)
#                                         ),
#                            strand
#                            )
#                    )
		#			sep="_")
