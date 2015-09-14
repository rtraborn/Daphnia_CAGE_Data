
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

Dp.bam <- c("/home/rtraborn/Daphnia/CAGE/TCO/Dp_CAGE_all_sorted.bam")

Dp.bai <- c("/home/rtraborn/Daphnia/CAGE/TCO/Dp_CAGE_all_sorted.bam.bai")

broad_gr <- GRanges("scaffold_1", IRanges(2092131,2095168),strand="+")

broad_gr

myFlag <- scanBamFlag()

Broad_param <- ScanBamParam(flag=myFlag,what=c("flag", "mrnm", "mpos"), which=GRanges("scaffold_1",IRanges(2092131,2095168),strand="+"))

aln <- readGAlignments(Dp.bam, index=Dp.bai,param=Broad_param,use.names=TRUE)

endDp <- start(aln)

strandDp <- strand(aln)

for (i in 1:length(endDp)) { endDp[i] + 1 -> endDp[i] }

myStrand <- GRanges("scaffold_1",IRanges(start(aln),endDp),strand=strandDp,names=names(aln))

Dp_align <- as(myStrand, "GAlignments")

is(Dp_align)

Dp_gr

gr_gene <- subsetByOverlaps(Dp_gr,broad_gr)

gr_gene

gr_gene2 <- gr_gene[gr_gene$name %in% "CDS"]

#gr_gene3 <- gr_gene[gr_gene$name %in% "five_prime_utr"]

#gr_gene4 <- append(gr_gene2, gr_gene3)

#bam <- scanBam(Dp.bam,index=Dp.bai)

gr_gene2

plot <- autoplot(Dp_align, which = broad_gr,method="raw",geom="line", color="#510099",stat="coverage")

is(plot)

p1 <- plot + coord_cartesian(xlim = c(2092131,2095168)) + theme_bw()

p2 <- ggplot() + geom_alignment(gr_gene2,type="exon") + coord_cartesian(xlim = c(2092131,2095168),ylim=c(0,1000)) + scale_x_continuous() + 
                             theme_bw()

labeled(p1) <- FALSE

labeled(p2) <- FALSE

tracks(p1,p2,heights=c(8,1)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
    panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave(file="broad_TSR.png", width=6, height=2,dpi=300,scale=1.5)

peaked_gr <- GRanges("scaffold_4", IRanges(2460951,2463321),strand="+")

gr_gene3 <- subsetByOverlaps(Dp_gr,peaked_gr)

gr_gene3

gr_gene4 <- gr_gene3[gr_gene3$name %in% "CDS"]

gr_gene4

peaked_gr

myFlag <- scanBamFlag()

peaked_param <- ScanBamParam(flag=myFlag,what=c("flag", "mrnm", "mpos"), which=GRanges("scaffold_4",IRanges(2460951,2462370),strand="+"))

peaked_param

Broad_param

p_aln <- readGAlignments(Dp.bam, index=Dp.bai,param=peaked_param,use.names=TRUE)

endDp <- start(p_aln)

strandDp <- strand(p_aln)

strandDp

for (i in 1:length(endDp)) { endDp[i] + 1 -> endDp[i] }

myStrand <- GRanges("scaffold_4",IRanges(start(p_aln),endDp),strand=strandDp,names=names(p_aln))

Dp_align_p <- as(myStrand, "GAlignments")

is(Dp_align_p)

plot_peaked <- autoplot(Dp_align_p,which=peaked_gr,method="raw",geom="line",color="#E5A118",stat="coverage")

p1_p <- plot_peaked + coord_cartesian(xlim = c(2460951,2463321),ylim=c(0,1000)) + theme_bw()

p2_p <- ggplot() + geom_alignment(gr_gene4) + coord_cartesian(xlim = c(2460951,2463321)) + scale_x_continuous() + 
                             theme_bw()

tracks(p1_p,p2_p,heights=c(8,1)) + theme(panel.grid.major = element_blank(),
                                   panel.grid.minor = element_blank(), 
    panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave(file="peaked_TSR.png", width=6, height=2,dpi=300,scale=1.5)
