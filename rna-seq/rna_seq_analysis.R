#Project: Romee Lab CD25 Project
#Author: Mohamad Najia
#Objective: analyze RNA-Seq on CD4/8 CAR-T cells isolated from tumor-bearing mice

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
#Function Declarations
####################################################

scaleRows <- function(x) {
  rm <- rowMeans(x)
  x <- sweep(x, 1, rm)
  sx <- apply(x, 1, sd)
  x <- sweep(x, 1, sx, "/")
}


####################################################
#Set-up Environment
####################################################

#declare colormaps
pal_rna <- colorRampPalette(c("#352A86","#343DAE","#0262E0","#1389D2",
                              "#2DB7A3","#A5BE6A","#F8BA43","#F6DA23","#F8FA0D"))(100)
pal_atac <- colorRampPalette(c('#3361A5', '#248AF3', '#14B3FF', 
                               '#88CEEF', '#C1D5DC', '#EAD397', 
                               '#FDB31A','#E42A2A', '#A31D1D'))(100)

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



####################################################
#DEG Analysis
####################################################

tcell_types <- unique(df.samples$cell_type)
y_max <- c(10,30)
names(y_max) <- tcell_types

for (i in tcell_types) {
  #filter samples
  df.subset <- filter(df.samples, cell_type == i)
  
  #import kallisto abundances
  tsv_files <- paste0(tsv_dir, df.subset$sample_name, "/KALLISTO/abundance.tsv")
  names(tsv_files) <- df.subset$sample_name
  txi <- tximport(tsv_files, type = "kallisto", tx2gene = g2t_map)
  
  #set up DESeq 
  sampleTable <- df.subset[,c("donor", "treatment"), drop = FALSE]
  sampleTable$treatment <- factor(sampleTable$treatment, levels = c("CD19-CAR-CD25", "CD19-CAR"))
  
  dds <- DESeqDataSetFromTximport(txi, sampleTable, design = ~ donor + treatment)
  dds <- DESeq(dds)
  
  #perform PCA
  nrow(dds)
  keep <- rowSums(counts(dds)) > 1
  dds.keep <- dds[keep,]
  nrow(dds.keep)
  
  #vsd <- vst(dds.keep, blind = FALSE)
  #head(assay(vsd), 3)
  #colData(vsd)
  #pcaData <- plotPCA(vsd, intgroup = c("treatment"), returnData = TRUE)
  
  rld <- rlog(dds.keep, blind = FALSE)
  pcaData <- plotPCA(rld, intgroup = c("treatment"), returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  pcaData$donor <- df.cd8$donor
  
  gg <- ggplot(pcaData, aes(x=PC1, y=PC2, fill=treatment, color=treatment)) + 
    geom_point(aes(shape = factor(donor))) +
    xlab(paste0("PC1: (", percentVar[1], "% Variance Explained)")) +
    ylab(paste0("PC2: (", percentVar[2], "% Variance Explained)")) +
    ggtitle(paste0(i, " cells RNA-seq PCA")) + 
    theme_bw() +
    theme(axis.text.x = element_text(colour = "black"),
          axis.text.y = element_text(colour = "black"),
          plot.title = element_text(hjust = 0.5),
          legend.title = element_blank())
  
  pdf(file = paste0(output_dir, i, "_pca_rlog_trandformed_RNA_seq_data.pdf"), width = 4.5, height = 3, useDingbats = FALSE)  
  print(gg)
  dev.off()
  
  #extract comparisons for differential expression analysis
  resSig <- data.frame(results(dds, contrast=c("treatment", "CD19-CAR-CD25", "CD19-CAR"), parallel = TRUE, pAdjustMethod = "fdr", alpha = 0.05))
  resSig %>% filter(complete.cases(.)) %>% arrange(-log2FoldChange) %>% filter(padj < 0.05) -> out
  out$padj_nlog10 <- -log10(out$padj)
  out$gene <- rownames(out)
  
  #round
  out$baseMean <- round(out$baseMean, 2)
  out$log2FoldChange <- round(out$log2FoldChange, 2)
  out$pvalue <-sprintf("%.3e", out$pvalue)
  out$padj <-sprintf("%.3e", out$padj)
  
  #output the differential expression table for statistically significant genes
  write.table(out, 
              file = paste0(output_dir, "DESeq/", i, "_CD19-CAR-CD25_v_CD19-CAR_Padj05.tsv"), 
              row.names = FALSE, 
              col.names = TRUE, 
              sep = "\t", 
              quote = FALSE)
  
  #output the differential expression table for all genes (significant or not)
  resSig %>% filter(complete.cases(.)) %>% arrange(-log2FoldChange) -> out.complete
  out.complete$padj_nlog10 <- -log10(out.complete$padj)
  out.complete$gene <- rownames(out.complete)
  
  write.table(out.complete, 
              file = paste0(output_dir, "DESeq/", i, "_CD19-CAR-CD25_v_CD19-CAR.tsv"), 
              row.names = FALSE, 
              col.names = TRUE, 
              sep = "\t", 
              quote = FALSE)
  
  #plot volcano of DEGs
  df.plot <- out.complete
  df.plot$significant <- FALSE
  df.plot[(df.plot$padj < 0.05), "significant"] <- TRUE
  df.plot$label <- ""
  df.plot[df.plot$significant, "label"] <- df.plot[df.plot$significant, "gene"]
  
  gg <- ggplot(df.plot, aes(x=log2FoldChange, y=padj_nlog10, color=significant, alpha=significant)) + 
    geom_point() + 
    scale_color_manual(values=c("#999999", "#E69F00")) + 
    scale_alpha_manual(values=c(0.2,1)) + 
    scale_size_manual(values=c(0.75,1)) + 
    geom_hline(yintercept=(-1*log10(0.05)), linetype="dashed", color = "black") + 
    geom_vline(xintercept=0, linetype="dashed", color = "black") + 
    scale_x_continuous(limits = c(-10, 10)) + 
    scale_y_continuous(limits = c(0, y_max[i])) + 
    xlab("log2(fold change)") + 
    ylab("-log10(P-adjusted)") + 
    ggtitle(paste0(i, " cells DEGs") ) + 
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(colour="black"),
          legend.position = "none",
          axis.text.y = element_text(colour="black"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    geom_text_repel(data = df.plot,
                    aes(label = label), 
                    segment.color = "black",
                    colour = "black",
                    size = 4,
                    box.padding = 1,
                    max.overlaps = 500,
                    show.legend = FALSE)
  
  pdf(file = paste0(output_dir, i, "_DEG_volcano_plot.pdf"), width = 6, height = 6, useDingbats = FALSE)
  print(gg)
  dev.off()
  
}



####################################################
#DEG Interpretation: GSEA using Hallmark gene sets
####################################################
library(msigdbr)
library(org.Hs.eg.db)
library(AnnotationDbi)

m_df <- msigdbr(species = "Homo sapiens")
m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene)

#CD8T cells 
fn <- paste0(output_dir, "DESeq/CD8T_CD19-CAR-CD25_v_CD19-CAR.tsv")
#fn <- paste0(output_dir, "DESeq/CD4T_CD19-CAR-CD25_v_CD19-CAR.tsv")
currset <- fread(fn, data.table = FALSE)

currset <- currset %>% mutate(ENTREZID = mapIds(org.Hs.eg.db, keys = gene, column = "ENTREZID", keytype = "SYMBOL") %>% unname())
geneList <- setNames(currset$log2FoldChange,as.numeric(unlist(currset$ENTREZID)) )

em <- GSEA(geneList, TERM2GENE = m_t2g)
df.gsea <- em@result
dotplot(em, showCategory=30)

pdf(paste0(output_dir, "CD8T_GSEA_Hallmark_MYC_Targets_v1_NES_-1.76_padj_0.00527.pdf"), width = 3.5, height = 3.5, useDingbats = FALSE)
gseaplot2(em, 
          geneSetID = 1, 
          title = "CD8 T CD19-CAR-CD25 v CD19-CAR", 
          color = "black", 
          pvalue_table = TRUE,
          rel_heights = c(1,0.1,.4), 
          base_size = 11)
dev.off()




####################################################
#DEG Interpretation: Pathway enrichment
####################################################

#pathway enrichment in CD8s
fn <- paste0(output_dir, "kallisto_rna_seq_TPM_matrix.tsv")
df <- fread(fn, data.table = FALSE)
rownames(df) <- df$V1

geneUniverse <- rownames(df$V1)
geneUniverse <- unlist(mget(geneUniverse, envir=org.Hs.egALIAS2EG,
                            ifnotfound = NA))

#import DEGs
fn <- paste0(output_dir, "DESeq/CD8T_CD19-CAR-CD25_v_CD19-CAR.tsv")
currset <- fread(fn, data.table = FALSE)

#get up-regulated genes
deGenes <- currset %>% filter(log2FoldChange > 0)
#deGenes <- deGenes$gene

#check that the genes are found within the pathway database 
deGenes <- unlist(mget(deGenes$gene, envir=org.Hs.egALIAS2EG,
                       ifnotfound = NA))

#perform enrichment
ans.go <- enrichGO(gene = deGenes, 
                   ont = "BP",#"BP", #"MF", "CC"
                   OrgDb ="org.Hs.eg.db",
                   universe = geneUniverse,
                   readable = TRUE,
                   pvalueCutoff = 0.05)

#plot the pathways
#pdf(paste0(pathway_dir, comparison, "_upregulated_genes_enriched_pathways.pdf"), width = 5, height = 5, useDingbats = FALSE)
dotplot(ans.go) + ggtitle("GO Biological Process")
#dev.off()

#filter the significant pathways and plot 
tab.go <- as.data.frame(ans.go)


gene_list <- paste0(tab.go$geneID, collapse = "/") %>% str_split(pattern = "/") 
unique_gene_list <- gene_list[[1]] %>% unique()

tab.go <- subset(tab.go, Count>5)
tab.go$logFDR <- -log10(tab.go$p.adjust)
tab.go <- tab.go[1:15,]
tab.go$Description <- factor(tab.go$Description, levels=tab.go$Description)

gg <- ggplot(tab.go, aes(y=logFDR, x=Description, size=logFDR, color=Count)) + 
  geom_point() + 
  #geom_bar(stat="identity", color="black", fill = "lightblue") +
  coord_flip() + 
  xlab("") + 
  ylab("-log10(FDR)") + 
  #ggtitle("CD8 T cells: CD19-CAR-CD25 Downregulated Genes") + 
  scale_x_discrete(limits = rev(levels(tab.go$Description))) + 
  labs(color='Gene Count', size = "-log10(FDR)") +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        #legend.title = element_blank(), 
        #legend.position = "bottom", 
        plot.title = element_text(hjust = 0.5))




#KEGG Pathways

#import DEGs
fn <- paste0(output_dir, "DESeq/CD8T_CD19-CAR-CD25_v_CD19-CAR.tsv")
currset <- fread(fn, data.table = FALSE)
rownames(currset) <- currset$gene

#check that the genes are found within the pathway database 
geneIDs <- unlist(mget(currset$gene, envir=org.Hs.egALIAS2EG,
                       ifnotfound = NA))

geneIDs <- geneIDs[!is.na(geneIDs)]

deGenes <- currset[names(geneIDs), "log2FoldChange"]
names(deGenes) <- geneIDs

deGenes <- sort(deGenes[!is.na(deGenes)], decreasing = TRUE)

kk2 <- gseKEGG(geneList     = deGenes,
               organism     = 'hsa',
               minGSSize    = 20,
               pvalueCutoff = 0.05,
               verbose      = FALSE)

dotplot(kk2) + ggtitle("KEGG Pathways")

df.gsea <- kk2@result

df.plot <- df.gsea[, c("Description", "NES")]
df.plot$direction <- "up"
df.plot[df.plot$NES < 0, "direction"] <- "down"

df.plot <- df.plot[order(df.plot$NES, decreasing = FALSE),]
df.plot <- df.plot[df.plot$Description %in% c("Cell adhesion molecules", "Ribosome", "Antigen processing and presentation", "HIF-1 signaling pathway"),]
df.plot$Description <- factor(df.plot$Description, levels = df.plot$Description)

pdf(file = paste0(output_dir, "CD8T_GSEA_KEGG_pathways_CD19-CAR_v_CD19-CAR-CD25.pdf"), width = 4, height = 1.5)  

ggplot(df.plot, aes(x = NES, y = Description, fill = direction)) + 
  geom_bar(stat = "identity") +
  scale_x_continuous(limits = c(-2.5, 2.5)) + 
  xlab("NES") + 
  ylab("KEGG Pathways\n(FDR < 0.05)") +
  ggtitle("CD8 T cells GSEA on KEGG Pathways") + 
  theme_bw() +
  theme(legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour="black"), 
        axis.text.y = element_text(colour="black"))

dev.off()




