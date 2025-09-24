#Project: Romee Lab CD25 Project
#Objective: infer TF activity with chromVAR using ATAC-seq data
#Author: Mohamad Najia

library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(irlba)
library(annotables)
library(BuenColors)
library(chromVAR)
library(ChrAccR)
library(SummarizedExperiment)
library(BiocParallel)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(circlize)
library(ComplexHeatmap)
library(GenomicRanges)
library(ggplot2)
library(ggrepel)
register(MulticoreParam(2))
set.seed(14651)


############ Functions ############ 
scaleRows <- function(x) {
  rm <- rowMeans(x)
  x <- sweep(x, 1, rm)
  sx <- apply(x, 1, sd)
  x <- sweep(x, 1, sx, "/")
}

calcDiffMotifScores <- function(sgT, sgNTC) {
  data.frame(diff_score = sgT - sgNTC) %>% arrange(desc(diff_score)) -> motifScoreDiff
  motifScoreDiff$rank <- 1:dim(motifScoreDiff)[1]
  
  anames <- str_split_fixed(rownames(motifScoreDiff), pattern = "\\|", n = 3)
  motifScoreDiff$name <- paste0("(", anames[,1], ") ", anames[,2])
  
  return(motifScoreDiff)
}
###################################



### Initialize and Import Data ###

#declare colormaps
pal_rna <- colorRampPalette(c("#352A86","#343DAE","#0262E0","#1389D2",
                              "#2DB7A3","#A5BE6A","#F8BA43","#F6DA23","#F8FA0D"))(100)
pal_atac <- colorRampPalette(c('#3361A5', '#248AF3', '#14B3FF', 
                               '#88CEEF', '#C1D5DC', '#EAD397', 
                               '#FDB31A','#E42A2A', '#A31D1D'))(100)

#set up environment
project_dir <- "/Volumes/broad_blainey_lab_storage/MN/20250122_AA_Tcell_in_vivo_SS2_ATAC_nextseq/CD25_CAR-T_Manuscript/atac-seq/"

#import count matrix
fn <- paste0(project_dir, "ATAC_seq_counts_matrix.tsv")
df.countMatrix <- fread(fn, data.table = FALSE)

df.peaks <- df.countMatrix[,c(1:3)]
df.countMatrix <- df.countMatrix[,c(4:9)]


### Apply chromVAR to all peaks ###

#filter peaks to include only standard chromosomes 
inds <- df.peaks$chr %in% paste0("chr", c(1:22, "X", "Y"))
gr.peaks <- makeGRangesFromDataFrame(df.peaks[inds,])
df.countMatrix <- df.countMatrix[inds,]

#sort the peaks
inds <- order(gr.peaks)
gr.peaks <- gr.peaks[inds]
df.countMatrix <- df.countMatrix[inds,]

#chromVAR setup
SE <- SummarizedExperiment(
  rowRanges =  gr.peaks,
  colData = data.frame(Sample = colnames(df.countMatrix), 
                       Construct = str_split_fixed(colnames(df.countMatrix), pattern = "_", n=4)[,4]), 
  assays = list(counts = as.matrix(df.countMatrix))
)

SE <- filterPeaks(SE)
SE <- addGCBias(SE, genome = BSgenome.Hsapiens.UCSC.hg38)

#match motifs to peaks using Vierstra motif archetype database
fn <- paste0(project_dir, "Vierstra_Archetype_Motifs_v2.1.rds")
varchetypes <- readRDS(fn)
motif_ix <- matchMotifs(varchetypes, SE, genome = BSgenome.Hsapiens.UCSC.hg38)

#calculate motif deviations
dev <- computeDeviations(object = SE, annotations = motif_ix)

#get the chromVAR deviations and z-scores
cvdevs <- deviations(dev)
cvzscores <- deviationScores(dev)

#determine differential deviations
difdev <- differentialDeviations(dev, "Construct")


#plot TF motif deviations for each donor

#Donor01
motifScoreDiff1 <- calcDiffMotifScores(cvzscores[,2], cvzscores[,1])
qbounds <- quantile(motifScoreDiff1$diff_score, c(0.25, 0.75))
donor1_bottom <- motifScoreDiff1[motifScoreDiff1$diff_score < qbounds[1], "name"]
donor1_top <- motifScoreDiff1[motifScoreDiff1$diff_score > qbounds[2], "name"]

#Donor02
motifScoreDiff2 <- calcDiffMotifScores(cvzscores[,4], cvzscores[,3])
qbounds <- quantile(motifScoreDiff2$diff_score, c(0.25, 0.75))
donor2_bottom <- motifScoreDiff2[motifScoreDiff2$diff_score < qbounds[1], "name"]
donor2_top <- motifScoreDiff2[motifScoreDiff2$diff_score > qbounds[2], "name"]

#Donor03
motifScoreDiff3 <- calcDiffMotifScores(cvzscores[,6], cvzscores[,5])
qbounds <- quantile(motifScoreDiff3$diff_score, c(0.25, 0.75))
donor3_bottom <- motifScoreDiff3[motifScoreDiff3$diff_score < qbounds[1], "name"]
donor3_top <- motifScoreDiff3[motifScoreDiff3$diff_score > qbounds[2], "name"]

#find common motifs across donors
common_bottom <- intersect(intersect(donor1_bottom, donor2_bottom), donor3_bottom)
common_bottom <- common_bottom[c(1,2,3,4,5,6,7,8,9,10,11,12)]
common_top <- intersect(intersect(donor1_top, donor2_top), donor3_top)
common_top <- common_top[c(1,2,3,4,5,9,10,11,12,13)]

motifScoreDiff1$label <- ""
motifScoreDiff1[motifScoreDiff1$name %in% c(common_bottom, common_top), "label"] <- motifScoreDiff1[motifScoreDiff1$name %in% c(common_bottom, common_top), "name"]

motifScoreDiff2$label <- ""
motifScoreDiff2[motifScoreDiff2$name %in% c(common_bottom, common_top), "label"] <- motifScoreDiff2[motifScoreDiff2$name %in% c(common_bottom, common_top), "name"]

motifScoreDiff3$label <- ""
motifScoreDiff3[motifScoreDiff3$name %in% c(common_bottom, common_top), "label"] <- motifScoreDiff3[motifScoreDiff3$name %in% c(common_bottom, common_top), "name"]


gg1 <- ggplot(motifScoreDiff1, aes(x = rank, y = diff_score, label = name)) + 
  geom_point(size = 0.75) +
  pretty_plot(fontsize = 8) + 
  L_border() + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
  xlab("Motif rank") + 
  ylab("Delta chromVAR deviation score") + 
  ggtitle("Donor 1") +
  theme_bw() + 
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() ) +
  geom_text_repel(data = motifScoreDiff1,
                  aes(label = label), 
                  segment.color = "black",
                  colour = "black",
                  size = 4,
                  box.padding = 1,
                  max.overlaps = 500,
                  show.legend = FALSE)
  #geom_label_repel(data = subset(motifScoreDiff, diff_score > 4),
  #                 box.padding   = 1, 
  #                 point.padding = 0.5,
  #                 segment.color = 'grey50',
  #                 max.overlaps = Inf) +
  #geom_label_repel(data = subset(motifScoreDiff, diff_score < -4.5),
  #                 box.padding   = 1, 
  #                 point.padding = 0.5,
  #                 segment.color = 'grey50',
  #                 max.overlaps = Inf)


gg2 <- ggplot(motifScoreDiff2, aes(x = rank, y = diff_score, label = name)) + 
  geom_point(size = 0.75) +
  pretty_plot(fontsize = 8) + 
  L_border() + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
  xlab("Motif rank") + 
  ylab("Delta chromVAR deviation score") + 
  ggtitle("Donor 2") +
  theme_bw() + 
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() ) +
  geom_text_repel(data = motifScoreDiff2,
                  aes(label = label), 
                  segment.color = "black",
                  colour = "black",
                  size = 4,
                  box.padding = 1,
                  max.overlaps = 500,
                  show.legend = FALSE)

gg3 <- ggplot(motifScoreDiff3, aes(x = rank, y = diff_score, label = name)) + 
  geom_point(size = 0.75) +
  pretty_plot(fontsize = 8) + 
  L_border() + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
  xlab("Motif rank") + 
  ylab("Delta chromVAR deviation score") + 
  ggtitle("Donor 3") +
  theme_bw() + 
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() ) +
  geom_text_repel(data = motifScoreDiff3,
                  aes(label = label), 
                  segment.color = "black",
                  colour = "black",
                  size = 4,
                  box.padding = 1,
                  max.overlaps = 500,
                  show.legend = FALSE)


#pdf(paste0(output_dir, "chromVAR_delta_dev_scores_vierstra_archetype_motifs_CD19-CAR-CD25_v_CD19-CAR.pdf"), width = 10, height = 6, useDingbats = FALSE)
gg1+gg2+gg3
#dev.off()


#take the average of delta chromVAR scores across all donors
donor1 <- motifScoreDiff1[order( motifScoreDiff1$name, decreasing = TRUE),]
donor2 <- motifScoreDiff2[order( motifScoreDiff2$name, decreasing = TRUE),]
donor3 <- motifScoreDiff3[order( motifScoreDiff3$name, decreasing = TRUE),]

tmp <- data.frame(donor1 = donor1$diff_score,
                  donor2 = donor2$diff_score,
                  donor3 = donor3$diff_score)


motifScoreDiff <- donor1
motifScoreDiff$diff_score <- rowMeans(tmp)

motifScoreDiff <- motifScoreDiff[order( motifScoreDiff$diff_score, decreasing = TRUE),]
motifScoreDiff$rank <- 1:dim(motifScoreDiff)[1]


gg <- ggplot(motifScoreDiff, aes(x = rank, y = diff_score, label = name)) + 
  geom_point(size = 0.75) +
  pretty_plot(fontsize = 8) + 
  L_border() + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
  xlab("Motif rank") + 
  ylab("Delta chromVAR deviation score") + 
  ggtitle("Differential TF motifs in CD8 T cells") +
  theme_bw() + 
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() ) +
  geom_text_repel(data = motifScoreDiff,
                  aes(label = label), 
                  segment.color = "black",
                  colour = "black",
                  size = 4,
                  box.padding = 1,
                  max.overlaps = 500,
                  show.legend = FALSE)

pdf(paste0(project_dir, "chromVAR_delta_dev_scores_vierstra_archetype_motifs_donor_avg_CD19-CAR-CD25_v_CD19-CAR.pdf"), width = 6, height = 6, useDingbats = FALSE)
gg
dev.off()





