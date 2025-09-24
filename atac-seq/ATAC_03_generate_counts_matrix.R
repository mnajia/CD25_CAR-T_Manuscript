#Project: Romee Lab CD25 Project
#Objective: generate a normalized ATAC-seq counts matrix 
#Author: Mohamad Najia

library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(GenomicRanges)
library(LSD)
library(preprocessCore)
library(edgeR)
library(ChIPseeker)



#declare variables
project_dir <- "/Volumes/broad_blainey_lab_storage/MN/20250122_AA_Tcell_in_vivo_SS2_ATAC_nextseq/CD25_CAR-T_Manuscript/atac-seq/"
counts_dir <- paste0(project_dir, "counts/")
consensus_peaks <- paste0(project_dir, "AA_PB_T_in_vivo_ATAC_hg38_consensusPeaks_minSamples2_minScore1.narrowPeak")
samplesheet <- paste0(project_dir, "samples_annotation.csv")

#import consensus peaks
df.peaks <- fread(consensus_peaks, header = FALSE, data.table = FALSE)
df.peaks <- df.peaks[,c(1:3)]
colnames(df.peaks) <- c("chr", "start", "end")

#import samplesheet 
samples <- fread(input = samplesheet, header = TRUE, data.table = FALSE)
countsList <- List()

#import individual sample counts 
for (sampleName in samples$sample_name) {
  print(sampleName)
  
  counts <- fread(input = paste0(counts_dir, sampleName, "_counts.bed"), header = FALSE, data.table = FALSE)
  countsList[[sampleName]] <- counts[,11]
}

#combine into a single matrix 
countsMatrix <- do.call("cbind", countsList)

#add peak information
df.counts <- cbind(df.peaks, countsMatrix)

#save the raw counts matrix 
write.table(df.counts,
            file = paste0(project_dir, "ATAC_seq_counts_matrix.tsv"),
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            sep = "\t"
            )

#perform CPM normalization 
cpm.df <- cpm(countsMatrix, log=TRUE, prior.count=2)

#Perform quantile normalization
norm.df <- normalize.quantiles(cpm.df)
norm.df <- as.data.frame(norm.df)
colnames(norm.df) <- colnames(cpm.df)

#save the normalized counts matrix 
df.norm <- cbind(df.peaks, norm.df)

write.table(df.norm,
            file = paste0(project_dir, "ATAC_seq_normalized_counts_matrix.tsv"),
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            sep = "\t"
)


#visualize replicate concordance as scatterplots 
heatscatter(norm.df$`ATAC_Donor01_CD8T_CD19-CAR`, norm.df$`ATAC_Donor02_CD8T_CD19-CAR`)
cor(norm.df$`ATAC_Donor01_CD8T_CD19-CAR`, norm.df$`ATAC_Donor02_CD8T_CD19-CAR`, method = "pearson")

pdf(file = paste(project_dir, "ATAC_pearson_cor.pdf", sep = "/"), width = 6, height = 5.5, useDingbats = FALSE)  
ComplexHeatmap::Heatmap( cor(norm.df, method = "pearson") )
dev.off()

library(MASS)
library(ggplot2)
library(viridis)

# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

dat <- norm.df[,c("ATAC_Donor01_CD8T_CD19-CAR","ATAC_Donor02_CD8T_CD19-CAR")]
colnames(dat) <- c("x","y")
dat$density <- get_density(dat$x, dat$y, n = 100)
ggplot(dat) + 
  geom_point(aes(x, y, color = density)) + 
  scale_color_viridis() + 
  theme_bw(base_size = 16)




