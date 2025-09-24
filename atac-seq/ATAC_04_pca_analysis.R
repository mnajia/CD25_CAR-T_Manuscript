#Project: Romee Lab CD25 Project
#Objective: perform PCA on a normalized ATAC-Seq counts matrix 
#Author: Mohamad Najia

library(irlba)
library(GenomicRanges)
library(matrixStats)
library(densityClust)
library(Rtsne)
library(ComplexHeatmap)
library(M3C)
library(ggplot2)
library(data.table)


#initialize
project_dir <- "/Volumes/broad_blainey_lab_storage/MN/20250122_AA_Tcell_in_vivo_SS2_ATAC_nextseq/CD25_CAR-T_Manuscript/atac-seq/"

#import samplesheet 
fn <- paste0(project_dir, "samples_annotation.csv")
samples <- fread(input = fn, data.table = FALSE)
rownames(samples) <- samples$sample_name

#import the counts matrix 
fn <- paste0(project_dir, "ATAC_seq_normalized_counts_matrix.tsv")
normMat <- fread(fn, data.table = FALSE)
peaks <- normMat[,c(1:3)]
normMat <- as.matrix(normMat[,c(4:9)])

#get top N peaks by variance across samples 
N <- 75000
df.topPeaks <- normMat[order(rowVars(normMat), decreasing = TRUE),][1:N,]

#perform PCA on top  matrix 
p1 <- prcomp_irlba(df.topPeaks, n = 5, scale. = FALSE)
samplePCs <- p1$rotation
cumulative_var <- summary(p1)$importance[3,]

#determine inflection point of PCs
df_plot <- data.frame(PC = 1:length(cumulative_var), var = cumulative_var)

#pdf(file = paste0(figures_dir, "normalized_count_matrix_20220301_cumulative_PCA_variance.pdf"), width = 4, height = 4, useDingbats = FALSE)  
ggplot(df_plot, aes(x=PC, y=var)) + geom_point() + 
  geom_vline(xintercept=9, linetype="dashed", color = "red") + 
  xlab("Principal Component") + 
  ylab("Cumulative Variance Explained") + 
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))
#dev.off()


df.pcs <- as.data.frame(samplePCs)
df.pcs$gRNA <- str_split_fixed(samples$sample_name, pattern = "_", n=4)[,4]
df.pcs$donor <- str_split_fixed(samples$sample_name, pattern = "_", n=4)[,2]

gg <- ggplot(df.pcs, aes(x=PC1, y=PC2, fill=gRNA, color=gRNA)) + 
  geom_point(aes(shape = factor(donor))) + 
  xlab("PC1 (76.2% Variance Explained)") + 
  ylab("PC2 (8.6% Variance Explained)") + 
  ggtitle("CD8 T cells ATAC-seq PCA") + 
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.title = element_blank())

pdf(file = paste(project_dir, "ATAC_pca_normalized_count_matrix.pdf", sep = "/"), width = 4.5, height = 3, useDingbats = FALSE)  
gg
dev.off()
