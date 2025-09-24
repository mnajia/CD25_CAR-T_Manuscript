#Project: Romee Lab CD25 Project
#Author: Mohamad Najia
#Objective: export kallisto RNA-seq expression matrices

library(DESeq2)
library(IHW)
library(jsonlite)
library(dplyr)
library(data.table)
library(tximport)
library(rhdf5)
library(ggplot2)
library(scales)
library(EnhancedVolcano)
library(rjson)
library(ComplexHeatmap)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(circlize)
library(irlba)
library(matrixStats)
library(M3C)
library(clusterProfiler)
library(enrichplot)



####################################################
#Set-up Environment
####################################################

#initialize variables
project_dir <- "/Volumes/broad_blainey_lab_storage/MN/20250122_AA_Tcell_in_vivo_SS2_ATAC_nextseq/CD25_CAR-T_Manuscript/rna-seq/"
output_dir <- paste0(project_dir, "rna_seq_analysis/")
tsv_dir <- paste0(project_dir, "kallisto_output/")

#import samplesheet
fn <- paste0(project_dir, "kallisto_samplesheet.txt")
df.samples <- fread(fn, header = FALSE, data.table = FALSE)
colnames(df.samples) <- c("sample_name", "fq1", "fq2")
rownames(df.samples) <- df.samples$sample_name

tmp <- str_split_fixed(df.samples$sample_name, pattern = "_", n=5)
df.samples$label <- paste0(tmp[,3], "_", tmp[,4])
df.samples$cell_type <- tmp[,3]
df.samples$donor <- tmp[,2]
df.samples$treatment <- tmp[,4]

#import hg19 USCS genomeStudio gene symbol to transcript id mappings
fn <- paste0(project_dir, "/kallisto_custom_index/hg19.annot.cdna.gene.symbol.transcript.map.CAR.custom")
g2t_map <- read.table(file = fn, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(g2t_map) <- c("transcript_id", "gene_symbol")

#import kallisto abundances
tsv_files <- paste0(tsv_dir, df.samples$sample_name, "/KALLISTO/abundance.tsv")
names(tsv_files) <- df.samples$sample_name
txi <- tximport(tsv_files, type = "kallisto", tx2gene = g2t_map)

#save kallisto expression matrices 
write.table(txi$abundance,
            file = paste0(output_dir, "kallisto_rna_seq_TPM_matrix.tsv"), 
            row.names = TRUE, 
            col.names = TRUE, 
            sep = "\t", 
            quote = FALSE)

write.table(txi$counts,
            file = paste0(output_dir, "kallisto_rna_seq_counts_matrix.tsv"), 
            row.names = TRUE, 
            col.names = TRUE, 
            sep = "\t", 
            quote = FALSE)




