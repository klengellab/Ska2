# The RNA-seq analysis of Ska2 knockdown in 2 and 4 weeks. 
# Authohr: Joy Otten Date: 12-01-2020

# Load in libraries
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(gplots)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(dplyr)
library(VennDiagram)
library(EnhancedVolcano)
library(topGO)
library(knitr)
library(kableExtra)
library(magrittr)
library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(sva)
library(biomaRt)
library(org.Mm.eg.db)
library(KEGGprofile)
options(stringsAsFactors = FALSE)
set.seed(1234)

# Functions:
ggPCA <- function(pca, metadata, variables){
  PCA_out <- as.data.frame(pca$x)
  percentage <- round(pca$sdev / sum(pca$sdev) * 100, 2)
  percentage <- paste0(colnames(PCA_out), " (", percentage, "%", ")")
  theme <- theme(panel.background = element_blank(), 
                 panel.border = element_rect(fill = NA), 
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(), 
                 strip.background = element_blank(), 
                 axis.text.x = element_text(colour = "black"), 
                 axis.text.y = element_text(colour = "black"),
                 axis.ticks = element_line(colour = "black"), 
                 plot.margin = unit(c(1, 1, 1, 1), "line"))
  for(i in variables){
    p <- ggplot(PCA_out, aes(x = PC1, y = PC2, color = metadata[, i])) + 
      geom_point() + theme + xlab(percentage[1]) + ylab(percentage[2]) + 
      labs(color = i)
    print(p)
  }
}
source("/PHShome/je637/gitlab/general_functions/plot/customPCA.R")
biomart_ensembl_to_gene <- function(gene){
  ensembl=useMart("ensembl")
  ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
  names <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                 filters = 'ensembl_gene_id', 
                 values = gene, 
                 mart = ensembl)
  return(names)
}

# Set Paths
output_path <- "/PHShome/je637/RNAseq/RNAseq_Ska2/output/"
data_path <- "/PHShome/je637/RNAseq/RNAseq_Ska2/data/"

# Load in data
data <- read.csv(paste0(data_path, "counts.csv"), header = TRUE)
metadata <- read.csv(paste0(data_path, "jakob_ska2.csv"), header = TRUE, sep = ",")

# The first column in the data are the ensemmbl gene names
rownames(data) <- data$gene
data <- data[,-1]

# Create metadata variabless
metadata$condition <- NA
x <- which(metadata$group == "SKA2-sh1" & metadata$timepoint == "2 weeks")
metadata[x,8] <- "SKA2_KD_week2"
x <- which(metadata$group == "SCR" & metadata$timepoint == "2 weeks")
metadata[x,8] <- "SCR_week2"
x <- which(metadata$group == "SKA2-sh1" & metadata$timepoint == "4 weeks")
metadata[x,8] <- "SKA2_KD_week4"
x <- which(metadata$group == "SCR" & metadata$timepoint == "4 weeks")
metadata[x,8] <- "SCR_week4"

metadata$num_condition <- as.numeric(as.factor(metadata$condition))
metadata$samplename <- paste0("X", metadata$samplename)
rownames(metadata) <- metadata$samplename

# Putting the metadata and expression data in the same order
x <- match(colnames(data),rownames(metadata))
metadata <- metadata[x,]
all(rownames(metadata) == colnames(data)) # This should be true if not check the order.

metadata$seq_depth <- colSums(data)
  
  
# Principal component analysis before normalization
PCA1 <- as.matrix(t(data))
PCA1 <- log2(PCA1 + 1)
PCA1 <- prcomp(PCA1)

ggPCA(PCA1, metadata, colnames(metadata)[3:10])
# The data groups based on the conditions

# custom PCA plot before normalization
pd <- metadata[,c(3:10)]
log_data <- log2(data + 1)
custom.PCA(beta = log_data, pd = pd, plot.title = "PCA before normalisation")
# There are significant co-variables related to sequencing depth. 

# Filtering counts that a gene needs a read with 10 or more counts
keep <- apply(data,1,function(x){any(x >=10)})
filtered_data <- data[keep, ]
filtered_data <- as.matrix(filtered_data)
dim(filtered_data) # 20819  32

barplot(colSums(filtered_data), main = "Total reads", xlab = "sample", 
        ylab = "log2 exprssionn filtered data")

log_data <- log2(filtered_data + 1)
custom.PCA(beta = log_data, pd = pd, plot.title = "PCA before normalisation")
# There are significant co-variables related to sequencing depth, lib batch, totalng and rnaconc.

# Differential expression analysis with Deseq2
coldata <- metadata[,c(1,7,8,9)]
coldata$libbatch <- as.factor(coldata$libbatch)
coldata$condition <- as.factor(coldata$condition)
all(rownames(coldata) == colnames(filtered_data)) # Make sure the output is true
dds <- DESeqDataSetFromMatrix(countData = filtered_data, colData = coldata, 
                              design = ~ condition + libbatch)
dds$condition <- factor(dds$condition, levels = c("SCR_week2", "SKA2_KD_week2", "SCR_week4", "SKA2_KD_week4"))dds <- DESeq(dds)
resultsNames(dds)

# Check the normalization
norm.data <-counts(dds, normalized=TRUE)
log_norm <- as.data.frame(log2(norm.data + 1))
barplot(colSums(log_norm), main = "Normalized Total reads", xlab = "sample", ylab = "log read counts")
pd <- metadata[,c(3:10)]
custom.PCA(beta = norm.data, pd = pd, plot.title = "PCA after normalization")
# Sequencing depth is not a significant co-variable anymore after normalization

# Results differential expression analysis week 2
Scr2vsSka2 <- results(dds, pAdjustMethod = "bonferroni", contrast=c("condition", "SKA2_KD_week2", "SCR_week2"))

sigScr2vsSka2 <- as.data.frame(Scr2vsSka2)
sigScr2vsSka2$genes <- rownames(Scr2vsSka2)
length(rownames(sigScr2vsSka2))
sig2 <- filter(sigScr2vsSka2, padj <= 0.05)
length(rownames(sig2))

## up and down regulated genes week 2
sigScr2vsSka2_up <- filter(sigScr2vsSka2, padj <= 0.05 & log2FoldChange > 2)
length(rownames(sigScr2vsSka2_up)) 
sigScr2vsSka2_down <- filter(sigScr2vsSka2, padj <= 0.05 & log2FoldChange < -2)
length(rownames(sigScr2vsSka2_down))

# Results differential expression analysis week 4
Scr4vsSka4 <- results(dds, pAdjustMethod = "bonferroni", contrast=c("condition", "SKA2_KD_week4", "SCR_week4"))
sigScr4vsSka4 <- as.data.frame(Scr4vsSka4)
sigScr4vsSka4$genes <- rownames(Scr4vsSka4)
length(rownames(sigScr4vsSka4)) 
sigScr4 <- filter(sigScr4vsSka4, padj <= 0.05)
length(rownames(sigScr4))
sigScr4vsSka4_up <- filter(sigScr4vsSka4, padj <= 0.05 & log2FoldChange > 2)
length(rownames(sigScr4vsSka4_up)) 
sigScr4vsSka4_down <- filter(sigScr4vsSka4, padj <= 0.05 & log2FoldChange < -2)
length(rownames(sigScr4vsSka4_down)) 

# To get external gene names from ensembl IDs
Ska2_2week <- as.list(as.vector(sig2$genes))
Ska2_4week <- as.list(as.vector(sigScr4$genes))

ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
x <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name'), 
           filters = 'ensembl_gene_id', 
           values = sig2$genes, 
           mart = ensembl)
y <- which(sig2$genes %in% x$ensembl_gene_id)
sig2 <- sig2[y,]
sig2$gene_id <- x$external_gene_name

x <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name'), 
           filters = 'ensembl_gene_id', 
           values = sigScr4$genes,
           mart = ensembl)
y <- which(sigScr4$genes %in% x$ensembl_gene_id)
sigScr4 <- sigScr4[y,]
sigScr4$gene_id <- x$external_gene_name

# Volcano plots of the differentially expressed genes
EnhancedVolcano(sigScr2vsSka2,
                lab = sigScr2vsSka2$genes,
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = NULL,
                xlim = c(-5, 8),
                title = 'Knockdown Ska2 2 weeks vs Control 2 weeks',
                subtitle = 'Differential expression',
                pCutoff = 10e-16,
                FCcutoff = 2,
                pointSize = 2.0,
                labSize = 0,
                colAlpha = 1)

EnhancedVolcano(sigScr4vsSka4,
                lab = sigScr4vsSka4$genes,
                x = 'log2FoldChange',
                y = 'pvalue',
                xlim = c(-5, 8),
                title = 'Knockdown Ska2 4 weeks vs Control 4 weeks',
                subtitle = 'Differential expression',
                pCutoff = 10e-16,
                FCcutoff = 2,
                pointSize = 2.0,
                labSize = 0,
                colAlpha = 1)


# Enrichment analysis

all_genes <- rownames(filtered_data)

x <- which(sigScr2vsSka2$padj < 0.05)
sig_week2 <- sigScr2vsSka2[x,]
x <- which(sigScr4vsSka4$padj < 0.05)
sig_week4 <- sigScr4vsSka4[x,]

## GO analysis 
enrich_2week <- enrichGO(
  gene = sig_week2$genes,
  OrgDb = "org.Mm.eg.db",
  keyType = "ENSEMBL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "bonferroni",
  universe = all_genes,
  qvalueCutoff = 0.05,
  readable = TRUE)
length(which(enrich_2week@result$p.adjust <= 0.05))
enrich_2week <- enrich_2week@result[c(1:342),]

enrich_4week <- enrichGO(
  gene = sig_week4$genes,
  OrgDb = "org.Mm.eg.db",
  keyType = "ENSEMBL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "bonferroni",
  universe = all_genes,
  qvalueCutoff = 0.05,
  readable = TRUE)
length(which(enrich_4week@result$p.adjust <= 0.05))
enrich_4week <- enrich_4week@result[c(1:153),]

# Kegg pathway analysis for the differential expressed genes in 2 and 4 weeks
search_kegg_organism('mmu', by='kegg_code')
mmus <- search_kegg_organism("Mus musculus", by="scientific_name")

# Convert Ensembl ID to NCBI ID
ensembl = useDataset("mmusculus_gene_ensembl", mart = ensembl)
# 2 weeks KEGG analysis
week2_KEGG <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                    filters = "ensembl_gene_id",
                    values = sig_week2$genes,
                    mart = ensembl)

# 4 weeks KEGG analysis
week4_KEGG <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                    filters = "ensembl_gene_id",
                    values = sig_week4$genes,
                    mart = ensembl)
all_genes_KEGG <- getBM(attributes=c("ensembl_gene_id", "entrezgene_id"),
                        filters = "ensembl_gene_id",
                        values = all_genes,
                        mart = ensembl)
all_genes_KEGG$entrezgene_id <- as.character(all_genes_KEGG$entrezgene_id)

# 2 weeks KEGG
kegg_2weeks <- enrichKEGG(gene = week2_KEGG$entrezgene_id, organism = "mmu", keyType = "kegg", 
                          pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = all_genes_KEGG$entrezgene_id)
kegg_2weeks <- as.data.frame(kegg_2weeks@result)

# 4 weeks KEGG
kegg_4weeks <- enrichKEGG(gene = week4_KEGG$entrezgene_id, organism = "mmu", keyType = "kegg", 
                          pvalueCutoff = 0.05, pAdjustMethod = "BH", universe = all_genes_KEGG$entrezgene_id)
kegg_4weeks <- as.data.frame(kegg_4weeks@result)


