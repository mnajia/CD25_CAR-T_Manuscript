#! /bin/bash
#manually generate pepatac consensus peakset
#Mohamad Najia


use UGER
ish -l h_vmem=20G

use Anaconda3
source activate /broad/blainey_lab/Mo/pepatac_env


R
devtools::install_github("databio/pepatac", subdir="PEPATACr")
library(PEPATACr)
pep <- "/broad/blainey_lab_storage/MN/20250122_AA_Tcell_in_vivo_SS2_ATAC_nextseq/project_config.yaml"

# Load the project metadata
prj <- pepr::Project(pep)
project_name    <- pepr::config(prj)$name
project_samples <- pepr::sampleTable(prj)$sample_name
sample_table    <- data.table::data.table(
    sample_name=pepr::sampleTable(prj)$sample_name,
    genome=pepr::sampleTable(prj)$genome)

# Specify file locations
output_dir  <- "/broad/blainey_lab_storage/MN/20250122_AA_Tcell_in_vivo_SS2_ATAC_nextseq/pepatac_pipeline"
results_dir <- file.path(output_dir, "results_pipeline")
summary_dir <- file.path(output_dir, "summary")

# Get project assets
missing_files <- 0
assets <- data.table(sample_name=character(), asset=character(), path=character(), annotation=character())

for (sample in project_samples) {
    sample_output_folder <- file.path(results_dir, sample)
    sample_assets_file   <- file.path(sample_output_folder, "assets.tsv")

    if (!file.exists(sample_assets_file)) {
        missing_files <- missing_files + 1
        next
    }

    t <- fread(sample_assets_file, header=FALSE, col.names=c('asset', 'path', 'annotation'))
    t <- t[!duplicated(t[, c('asset', 'path', 'annotation')], fromLast=TRUE),]
    t[,sample_name:=sample]
    assets = rbind(assets, t)
}

# Generate consensus peaks and write to project output directory
min_samples <- 1
min_score <- 1
peak_filepath <- consensusPeaks(sample_table, summary_dir, results_dir, assets, min_samples, min_score)

# Load the peak file into R
peaks <- data.table::fread(peak_filepath[[1]])


