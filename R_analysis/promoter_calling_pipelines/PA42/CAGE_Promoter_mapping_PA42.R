#A workflow to caculate the positions of TSRs (putative promoters) within our Daphnia pulex CAGE datasets using the CAGEr package. 
#Based on the vignette of the 'CAGEr' Bioconductor Package (Haberle et al., 2014)
#Adapted for Algning against the new PA42 assembly

require(CAGEr)
#load the genome as appropriate
require(BSgenome.Dpulex.IU.dappa2)
#enter the location of the directory of mapped bam files from CAGE or other 5' end reads. (Note: the directory must ONLY have the bam files of interest in it)
thisDir <- c("/home/rtraborn/Daphnia/CAGE/demultiplexed_matched/bam_files/PA42_bams/filtered/replicate_bams")
pathsToInputFiles <- list.files(thisDir, full.names = TRUE)   

#creates the CAGEset S4 object. Enter the genomeName and sampleLabels as appropriate for your analysis
myCAGEset <- new("CAGEset",genomeName="BSgenome.Dpulex.IU.dappa2",inputFiles=pathsToInputFiles,inputFilesType="bam",sampleLabels=c("mat_females_1","mat_females_2","mat_females_3","mat_males_1","mat_males_2","pE_females_1","pE_females_2","pE_females_3"))

#calculates and plots the pearson correlation between all samples
#corr.m <- plotCorrelation(myCAGEset, samples = "all", method = "pearson")

#calculating the TSS clusters (a subunit of TSRs)
getCTSS(myCAGEset)
ctss <- CTSStagCount(myCAGEset)
nCTSS <- nrow(ctss)
#outputs number of CTSS
cat("The number of CTSS in sample is", nCTSS,"\n")

#plots the reverse cumulative distribution of the CTSS in the sample
#in our case the alpha is set to 1.05 because we already know the values are around 1
plotReverseCumulatives(myCAGEset, fitInRange = c(5, 1000), onePlot = TRUE)
normalizeTagCount(myCAGEset, method = "simpleTpm")
#uncomment the following line if you want a bedgraph file of the CTSSs in the sample
#exportCTSStoBedGraph(myCAGEset, values = "normalized", oneFile = TRUE)

#saves an R binary of the data. 
save.image("Dp_PA42.RData")

#clustering the CTSS into so-called Tag Clusters/TCs, which we call TSRs
clusterCTSS(object = myCAGEset, threshold = 1, thresholdIsTpm = TRUE,nrPassThreshold = 1, method = "distclu", maxDist = 20,removeSingletons = TRUE, keepSingletonsAbove = 5, useMulticore = T, nrCores = 6)

##TSR widths and summary statistics
cumulativeCTSSdistribution(myCAGEset, clusters = "tagClusters")
quantilePositions(myCAGEset, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)

#aggregates tag clusters and then calculates the consensus clusters (promoters across all three conditions) within the sample
aggregateTagClusters(myCAGEset, tpmThreshold = 5, qLow = 0.1, qUp = 0.9, maxDist = 100)
cumulativeCTSSdistribution(myCAGEset, clusters = "consensusClusters")
quantilePositions(myCAGEset, clusters="consensusClusters",qLow=0.1,qUp=0.9)

consCl <- consensusClusters(myCAGEset)
write.table(consCl, file="PA42_consensus_clusters_tpm.txt", col.names=TRUE, row.names=FALSE, quote=FALSE)

TSR_summary <- summary(consCl)
write.table(TSR_summary,file="TSR_interquantile_summary_PA42.txt",sep=" ",col.names=TRUE,row.names=FALSE,quote=FALSE)

#exports a bed file of the TSRs' interquantile widths
exportToBed(myCAGEset, what = "consensusClusters", oneFile = TRUE)

#saves the CAGEset to a binary 'RData' file in your working directory
save(myCAGEset,file="Dp_CAGEset_PA42.RData")

#analysis is complete
print("Analysis is complete!")
q()
