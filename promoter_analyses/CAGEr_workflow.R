# Analysis of Daphnia CAGE data, including promoter annotation using CAGEr and differential expression analysis
# Credit to Dave Tang for his github tutorial (located here: https://github.com/davetang/ccl2/blob/master/analysis.R)
#Written on January 20, 2015

#Session Info

sessionInfo()

#Package installation
required_packages <- c('CAGEr','BSgenome.Dpulex','biomaRt','edgeR','GenomicFeatures','GO.db','GOstats','Sushi','multicore')

#install check
package_check <- required_packages %in% installed.packages()[,"Package"]

if(all(package_check[1:9])==FALSE){
    source("http://bioconductor.org/biocLite.R")
}

#quick package installation function
install_package <- function(x){
    biocLite(x)
}

#this installs the missing Bioconductor packages
sapply(required_packages[which(!package_check[1:9])], install_package)

#if it's missing, install the multicore package
if(!package_check[9]){
    install.packages("multicore",repos='http://cran.us.r-project.org')
}

#loads Bioconductor
#source("http://bioconductor.org/biocLite.R")
#biocLite()

#loads CAGEr
library(CAGEr)

#loads D. pulex genome BS object
library(BSgenome.Dpulex.JGI.dpulex)

#location of bamfiles (move files into this new location on your system)
inputDir <- "/DATA/GROUP/rtraborn/Daphnia_Project/Dp_data/merged_bamfiles"

pathsToInputFiles <- list.files(inputDir, full.names = TRUE)

#for debugging
pathsToInputFiles

#now, creating the CAGEr object

DpCAGE <- new("CAGEset",
              genomeName = "BSgenome.Dpulex.JGI.dpulex",
              inputFiles = pathsToInputFiles,
              inputFilesType = 'bam',
              sampleLabels = c('mat_females')
              )

#checking whether the DpCAGE object was successfully created
DpCAGE

getCTSS(DpCAGE, mappingQualityThreshold=10,sequencingQualityThreshold=-1)
DpCAGE_ctss <- CTSStagCount(DpCAGE)

library(multicore)

clusterCTSS(object = DpCAGE,
                        threshold = 1,
                        thresholdIsTpm = TRUE,
                        nrPassThreshold = 1,
                        method = "distclu",
                        maxDist = 20,
                        removeSingletons = TRUE,
                        keepSingletonsAbove = 5,
                        useMulticore = T,
                        nrCores = 6
                       )

#The cumulative sum of CAGE tags in the genome
              cumulativeCTSSdistribution(DpCAGE, clusters = "tagClusters", useMulticore = T, nrCores = 8)
              quantilePositions(DpCAGE, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)

              #aggregating tag clusters across multiple CAGE datasets
              aggregateTagClusters(DpCAGE, tpmThreshold = 5, qLow = 0.1, qUp = 0.9, maxDist = 100)

              #Extracting consensus clusters from CAGEset object
              Dp_CAGE_consensus <- consensusClusters(DpCAGE)

              #CAGE data based expression clustering
              getExpressionProfiles(DpCAGE, what = "Dp_CAGE_consensus",
                                                          tpmThreshold = 10, nrPassThreshold = 1,
                                                          method = "som", xDim = 5, yDim = 5
                                                          )

#let's see what this gives us:

##Reading in file: /home/rtraborn/R/x86_64-redhat-linux-gnu-library/3.1/CAGEr/merged_bamfiles/mat_fem_filtered_merged.bam...
#-> Filtering out low quality reads...
#Error in .Call(.NAME, ..., PACKAGE = PACKAGE) :
#      negative length vectors are not allowed
#Calls: getCTSS ... .XStringQualityToIntegerMatrix -> matrix -> unlist -> unlist -> .Call2 -> .Call
#Execution halted
#Potential solution found here: https://support.bioconductor.org/p/52717/
