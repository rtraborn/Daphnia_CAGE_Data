
library(ggbio)

library(rtracklayer)

library(GenomicAlignments)

library(Rsamtools)

library(biovizBase)

library(GenomicRanges)

Dp_genes <- "/DATA/GROUP/rtraborn/Daphnia_Project/Dp_data/FrozenGeneCatalog.bed"

Dp_fg <- read.table(file="/home/rtraborn/Daphnia/annotation_files/FrozenGeneCatalog20110204_unsorted.bed",header=FALSE)

Dp_fg <- Dp_fg[,-9:-10]

names(Dp_fg) <- c("chr","start","end","name2","score","strand","score2","name")

head(Dp_fg)

Dp_gr <- with(Dp_fg, GRanges(chr, IRanges(start, end), strand=strand,name=name))

Dp_gr

Dp.bam <- c("/home/rtraborn/Daphnia/CAGE/TCO/pE_fem_filtered_merged.bam")

Dp.bai <- c("/home/rtraborn/Daphnia/CAGE/TCO/pE_fem_filtered_merged.bam.bai")

coord_region <- GRanges("scaffold_47", IRanges(116142,127656),strand="-")

myFlag <- scanBamFlag()

#scaffold_47:151,539-157,785
my_param <- ScanBamParam(flag=myFlag,what=c("flag", "mrnm", "mpos"), which=GRanges("scaffold_47", IRanges(116142,127656),strand="-"))

aln <- readGAlignments(Dp.bam, index=Dp.bai,param=my_param,use.names=TRUE)

endDp <- end(aln)

strandDp <- strand(aln)

for (i in 1:length(endDp)) { endDp[i] - 1 -> endDp[i] }

myStrand <- GRanges("scaffold_47",IRanges(endDp,end(aln)),strand=strandDp,names=names(aln))

Dp_align <- as(myStrand, "GAlignments")

is(Dp_align)

Dp_gr

gr_gene <- subsetByOverlaps(Dp_gr,coord_region)

gr_gene

gr_gene2 <- gr_gene[gr_gene$name %in% "CDS"]

#gr_gene3 <- gr_gene[gr_gene$name %in% "five_prime_utr"]

#gr_gene4 <- append(gr_gene2, gr_gene3)

#bam <- scanBam(Dp.bam,index=Dp.bai)

gr_gene2

plot <- autoplot(Dp_align, which = coord_region, method="raw",geom="area", color="#CC79A7", coverage.fill="#CC79A7", stat="coverage")

is(plot)

p1 <- plot + coord_cartesian(xlim = c(116142,127656)) + theme_bw() + ylim(0,15000)

p2 <- ggplot() + geom_alignment(gr_gene2,type="exon") + coord_cartesian(xlim = c(116142,127656)) + scale_x_continuous() + theme_bw()

labeled(p1) <- FALSE

labeled(p2) <- FALSE

tracks(p1,p2,heights=c(8,1)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave(file="pE_region_VTG.png", width=6, height=2,dpi=300,scale=1.5)

######## mature females plot ###############

coord_region <- GRanges("scaffold_47", IRanges(116142,127656),strand="-")

gr_gene3 <- subsetByOverlaps(Dp_gr,coord_region)

gr_gene3

gr_gene4 <- gr_gene3[gr_gene3$name %in% "CDS"]

gr_gene4

peaked_gr <- gr_gene4

myFlag <- scanBamFlag()

matfem_fem_param<- ScanBamParam(flag=myFlag,what=c("flag", "mrnm", "mpos"), which=GRanges("scaffold_47",IRanges(116142,127656),strand="-"))

Dp.bam <- c("/home/rtraborn/Daphnia/CAGE/TCO/mat_fem_filtered_merged.bam")

Dp.bai <- c("/home/rtraborn/Daphnia/CAGE/TCO/mat_fem_filtered_merged.bam.bai")

p_aln <- readGAlignments(Dp.bam, index=Dp.bai,param=matfem_fem_param,use.names=TRUE)

endDp <- end(p_aln)

strandDp <- strand(p_aln)

strandDp

for (i in 1:length(endDp)) { endDp[i] - 1 -> endDp[i] }

myStrand <- GRanges("scaffold_47",IRanges(endDp,end(p_aln)),strand=strandDp,names=names(p_aln))

Dp_align_p <- as(myStrand, "GAlignments")

is(Dp_align_p)

plot_peaked <- autoplot(Dp_align_p,which=coord_region,method="raw",geom="area",color="darkgreen",fill="darkgreen",stat="coverage", ymin=0, ymax=20000 )

p1_p <- plot_peaked + coord_cartesian(xlim = c(116142,127656)) + theme_bw() + ylim(0,15000)

p2_p <- ggplot() + geom_alignment(gr_gene4) + coord_cartesian(xlim=c(116142,127656)) + scale_x_continuous() + 
                             theme_bw()

tracks(p1_p,p2_p,heights=c(8,1)) + theme(panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), 
    panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave(file="mat_fem_region_VTG.png", width=6, height=2,dpi=300,scale=1.5)

########## mature males plot ###########

coord_region <- GRanges("scaffold_47", IRanges(116142,127656),strand="-")

gr_gene3 <- subsetByOverlaps(Dp_gr,coord_region)

gr_gene3

gr_gene4 <- gr_gene3[gr_gene3$name %in% "CDS"]

gr_gene4

peaked_gr <- gr_gene4

myFlag <- scanBamFlag()

matmale_param <- ScanBamParam(flag=myFlag,what=c("flag", "mrnm", "mpos"), which=GRanges("scaffold_47",IRanges(116142,127656),strand="-"))

Dp.bam <- c("/home/rtraborn/Daphnia/CAGE/TCO/mat_male_filtered_merged.bam")

Dp.bai <- c("/home/rtraborn/Daphnia/CAGE/TCO/mat_male_filtered_merged.bam.bai")

p_aln <- readGAlignments(Dp.bam, index=Dp.bai,param=matmale_param,use.names=TRUE)

endDp <- start(p_aln)

strandDp <- strand(p_aln)

strandDp

for (i in 1:length(endDp)) { endDp[i] + 1 -> endDp[i] }

myStrand <- GRanges("scaffold_47",IRanges(endDp,end(p_aln)),strand=strandDp,names=names(p_aln))

Dp_align_p <- as(myStrand, "GAlignments")

is(Dp_align_p)

plot_peaked <- autoplot(Dp_align_p,which=coord_region,method="raw",geom="area",color="purple",fill="purple",stat="coverage")

p3_p <- plot_peaked + coord_cartesian(xlim = c(116142,127656)) + theme_bw() + ylim(0,15000)

p4_p <- ggplot() + geom_alignment(gr_gene4) + coord_cartesian(xlim=c(116142,127656)) + scale_x_continuous() + 
                             theme_bw()

tracks(p3_p,p4_p,heights=c(8,1)) + theme(panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), 
    panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave(file="mat_males_region_VTG.png", width=6, height=2,dpi=300,scale=1.5)

tracks(p1,p1_p,p3_p,p4_p,heights=c(5,5,5,2)) + theme(panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), 
    panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave(file="VTG_plot_combined.png", width=6, height=2, scale=1.5)


