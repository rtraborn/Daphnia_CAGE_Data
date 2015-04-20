#A workflow to caculate the positions of TSRs (putative promoters) within our Daphnia pulex CAGE datasets using the CAGEr package. 
#Based on the vignette of the 'CAGEr' Bioconductor Package (Haberle et al., 2014)

require(CAGEr)
#load the genome as appropriate
require(BSgenome.Dpulex.JGI.dpulex)
#enter the location of the directory of mapped bam files from CAGE or other 5' end reads. (Note: the directory must ONLY have the bam files of interest in it)
c("/home/rtraborn/Daphnia/CAGE/corrected_bamfile") -> thisDir
pathsToInputFiles <- list.files(thisDir, full.names = TRUE)   

#creates the CAGEset S4 object. Enter the genomeName and sampleLabels as appropriate for your analysis
myCAGEset <- new("CAGEset",genomeName="BSgenome.Dpulex.JGI.dpulex",inputFiles=pathsToInputFiles,inputFilesType="bam",sampleLabels=c("mat_fem_1","mat_fem_2","mat_fem_3","mat_males_1","mat_males_2","pE_females_1","pE_females_2","pE_females_3"))

#calculating the TSS clusters (a subunit of TSRs)
getCTSS(myCAGEset)
ctss <- CTSStagCount(myCAGEset)
nCTSS <- nrow(ctss)
#outputs number of CTSS
cat("The number of CTSS in sample is", nCTSS,"\n")

#merging the replicate CAGE libraries into samples from the same condition                 
mergeSamples(myCAGEset, mergeIndex = c(1,1,1,3,3,2,2,2), mergedSampleLabels = c("asexual_females","pE_females","sexual_males"))
                 
#plots the reverse cumulative distribution of the CTSS in your sample. You will have to go back and  adjust the alpha value on the plot that is generated to normalize your sample to the Power Law distribution. At present we have T equal to one million, so the program caluclates the expression values in 'tags per million (tpm)'
#in our case the alpha is set to 1.05 because we already know the values are around 1
plotReverseCumulatives(myCAGEset, fitInRange = c(5, 1000), onePlot = TRUE)
normalizeTagCount(myCAGEset, method = "powerLaw", fitInRange = c(5, 1000), alpha = 1.05, T = 1*10^6)
save(myCAGEset,file="Dp_CAGEset_pipe.RData")
#uncomment the following line if you want a bedgraph file of the CTSSs in the sample
#exportCTSStoBedGraph(myCAGEset, values = "normalized", oneFile = TRUE)

#clustering the CTSS into so-called Tag Clusters/TCs, which we call TSRs
clusterCTSS(object = myCAGEset, threshold = 1, thresholdIsTpm = TRUE,nrPassThreshold = 1, method = "distclu", maxDist = 20,removeSingletons = TRUE, keepSingletonsAbove = 5)
             
##TSR widths and summary statistics
cumulativeCTSSdistribution(myCAGEset, clusters = "tagClusters")
quantilePositions(myCAGEset, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)

#creates objects containing TCs with interquantile widths and writes the summary statistics to a table
asex_TC <- tagClusters(myCAGEset, sample = "asexual_females", returnInterquantileWidth = TRUE, qLow = 0.1, qUp = 0.9)
pE_TC <- tagClusters(myCAGEset, sample = "pE_females", returnInterquantileWidth = TRUE, qLow = 0.1, qUp = 0.9)
male_TC <- tagClusters(myCAGEset, sample = "sexual_males", returnInterquantileWidth = TRUE, qLow = 0.1, qUp = 0.9)             
TSR_summary1 <- summary(asex_TC)
TSR_summary2 <- summary(pE_TC)             
TSR_summary3 <- summary(male_TC)
write.table(TSR_summary1,file="asex_fem_TSR_interquantile_summary.txt",sep=" ",col.names=TRUE,row.names=FALSE)
write.table(TSR_summary2,file="pE_fem_TSR_interquantile_summary.txt",sep=" ",col.names=TRUE,row.names=FALSE)
write.table(TSR_summary3,file="males_TSR_interquantile_summary.txt",sep=" ",col.names=TRUE,row.names=FALSE)             

#exports a bed file of the TSRs' interquantile widths
exportToBed(object = myCAGEset, what = "tagClusters", qLow = 0.1, qUp = 0.9, oneFile = TRUE)

#saves the CAGEset to a binary 'RData' file in your working directory
save(myCAGEset,file="Dp_CAGEset_pipe.RData")

#analysis is complete
print("Analysis is complete!")
q()
