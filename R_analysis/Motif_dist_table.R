setwd("/home/rtraborn/Daphnia/Daphnia_CAGE_Data/promoter_analyses/de_novo_discovery/coverage_histograms")

library(ggplot2)

library(reshape2)

list.files()

DpTCO_dist <- read.table(file="DpTCO_dist_50_50_histogram_100_1bp_update",sep="\t",header=TRUE) 

DpTCO_dist2 <- read.table(file="DpTCO_dist_50_50_histogram_200_5bp",sep="\t",header=TRUE)

head(DpTCO_dist)

hist_loc <- DpTCO_dist2[,1]

#hist_loc2 <- DpTCO_dist2[,1]

motif1_dens <- DpTCO_dist2[,3]

motif2_dens <- DpTCO_dist2[,6]

motif3_dens <- DpTCO_dist2[,9]

motif4_dens <- DpTCO_dist2[,12]

motif5_dens <- DpTCO_dist2[,15]

motif6_dens <- DpTCO_dist2[,21]

motif7_dens <- DpTCO_dist2[,24]

options(scipen=6)

motif_hist <- cbind(motif4_dens,motif5_dens)

motif_hist_i <- cbind(motif1_dens,motif6_dens)

rownames(motif_hist) <- hist_loc

rownames(motif_hist_i) <- hist_loc

head(motif_hist)

head(motif_hist_i)

motif_hist_m <- melt(motif_hist)

motif_hist_i_m <- melt(motif_hist_i)

names(motif_hist_m) <- c("position","motif","density")

names(motif_hist_i_m) <- c("position","motif","density")

head(motif_hist_m)

p <- ggplot(motif_hist_m, aes(x=position, y=density),group=motif)

p + geom_line(aes(colour=motif)) + xlim(-100,100) + theme_bw() + labs(x="Position relative to TSS\n(in bp)",y="Density of sites")

ggsave(file="Dp_motif_density_TC_50_50_1bp_hist_100_100_5_Dpm_45.png")

p <- ggplot(motif_hist_i_m, aes(x=position, y=density),group=motif)

p + geom_line(aes(colour=motif)) + xlim(-100,100) + theme_bw() + labs(x="Position relative to TSS\n(in bp)",y="Density of sites")

ggsave(file="Dp_motif_density_TC_50_50_5bp_hist_100_100_Dpm_5_16.png")

hist_loc2 <- DpTCO_dist[,1]

motif1_dens <- DpTCO_dist[,3]

motif2_dens <- DpTCO_dist[,6]

motif3_dens <- DpTCO_dist[,9]

motif7_dens <- DpTCO_dist[,24]

motif_hist2 <- cbind(motif2_dens,motif3_dens)

motif_hist2 <- motif_hist2[-101,] #removing the zero row

motif_hist2

rownames(motif_hist2) <- hist_loc2[-101]

which(hist_loc2==0)

motif_hist_m2 <- melt(motif_hist2)

names(motif_hist_m2) <- c("position","motif","density")

p2 <- ggplot(motif_hist_m2, aes(x=position, y=density),group=motif)

p2 + geom_line(aes(colour=motif)) + xlim(-50,50) + theme_bw() + labs(x="Position relative to TSS\n(in bp)",y="Density of sites")

ggsave(file="Dp_motif_density_TC_50_50_1bp_hist_high_freq_zoom_Dpm2_3.png")

