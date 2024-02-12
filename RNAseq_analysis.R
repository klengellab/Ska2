# Ska2 2 and 4 week differential expression analysis 

# libraries
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
library(corrplot)
library(factoextra)
library(openxlsx)
library(msigdbr)
library(GOplot)
options(stringsAsFactors = FALSE)
set.seed(1234)

# Functions
source("/PHShome/je637/gitlab/rna-seq/Ska2_Jakob_2020/functions.R")
genes <- read.csv("/PHShome/je637/general/tables/ensembl_w_description.mouse.csv")

# Paths
output_path <- "/PHShome/je637/RNAseq/RNAseq_Ska2/output/" # set this to a path were you want the images/files stored
data_path <- "/PHShome/je637/RNAseq/RNAseq_Ska2/data/" # set this to the path were the data is stored

# Data
#proportions_celltype <- read.csv("/PHShome/je637/ska2/output/proportions_celltypes.csv")
metadata <- read.csv(paste0(output_path, "metadata.csv"))
all_cell_prop <- read.csv(paste0(output_path, "all_cell_type_prop.csv"))
data <- read.csv(paste0(data_path, "counts.csv"), header = TRUE)

# Filtering and PCA analysis of the data
# --------
# The first column in the data are the ensembl gene names
rownames(data) <- data$gene
data <- data[,-1]

# Putting the metadata and expression data in the same order
metadata <- metadata[,-1]
rownames(metadata) <- metadata$samplename
x <- match(colnames(data),rownames(metadata))
metadata <- metadata[x,]
all(rownames(metadata) == colnames(data)) # This should be true if not check the order.

# Add sequencing depth to metadata
metadata$seq_depth <- colSums(data)

# First I need to add the cell type proportions to the metadata therefore, I need to 
# first split the dataframe by cell types
cell_types <- unique(all_cell_prop$X2)
for(cell in cell_types){
  message(cell)
  x <- which(cell == all_cell_prop$X2)
  cells <- all_cell_prop[x,]
  reorder <- match(rownames(metadata), cells$X1)
  cells <- cells[reorder, ]
  if(cell == "NEURON"){
    metadata$prop_neuron <- NA
    metadata$prop_neuron <- cells$value
  }
  if(cell == "OLIGODENDROCYTE"){
    metadata$prop_oligo <- NA
    metadata$prop_oligo <- cells$value
  }
  if(cell == "ASTROCYTE"){
    metadata$prop_astrocyte <- NA
    metadata$prop_astrocyte <- cells$value
  }
  if(cell == "POLYDENDROCYTE"){
    metadata$prop_polydendrocyte <- NA
    metadata$prop_polydendrocyte <- cells$value
  }
  if(cell == "FIBROBLAST"){
    metadata$prop_fibroblast <- NA
    metadata$prop_fibroblast <- cells$value
  }
  if(cell == "MICROGLIA"){
    metadata$prop_microglia <- NA
    metadata$prop_microglia <- cells$value
  }
  if(cell == "ENDOTHELIAL"){
    metadata$prop_endothelial <- NA
    metadata$prop_endothelial <- cells$value
  }
  if(cell == "NEUROGENESIS"){
    metadata$prop_neurogenesis <- NA
    metadata$prop_neurogenesis <- cells$value
  }
  if(cell == "MURAL"){
    metadata$prop_mural <- NA
    metadata$prop_mural <- cells$value
  }
  if(cell == "EPENDYMA"){
    metadata$prop_ependyma <- NA
    metadata$prop_ependyma <- cells$value
  }
}

# Filtering counts that a gene needs a read with 10 or more counts
keep <- apply(data,1,function(x){any(x >=10)})
filtered_data <- data[keep, ]
filtered_data <- as.matrix(filtered_data)
dim(filtered_data)

# Perform PCA analysis with the cell types to see which ones are significant in co-variables
# and if there is an co-linearity
pd <- metadata[,c(3:20)]
log_data <- log2(filtered_data + 1)
custom.PCA(beta = log_data, pd = pd, plot.title = "PCA before normalization")

# --------

# Differential expression analysis with Deseq2

# --------
# Differential expression analysis with Deseq2 including the cell type proportions
coldata <- metadata[,c(1,3:20)]
coldata$libbatch <- as.factor(coldata$libbatch)
rownames(coldata) <- coldata$samplename

# The neuronal and microglia population are co-linear with one another

## including microglia population
dds_cell_micro <- DESeqDataSetFromMatrix(countData = filtered_data, colData = coldata,
                                       design = ~ libbatch + prop_microglia + condition)
dds_cell_micro$condition <- factor(dds_cell_micro$condition, levels = c("SCR_week2", "SKA2_KD_week2", "SCR_week4", "SKA2_KD_week4"))
dds_cell_micro <- DESeq(dds_cell_micro)
resultsNames(dds_cell_micro)

Scr2vsSka2_micro <- results(dds_cell_micro, pAdjustMethod = "bonferroni", contrast=c("condition", "SKA2_KD_week2", "SCR_week2"))
Scr2vsSka2_micro <- as.data.frame(Scr2vsSka2_micro)
Scr2vsSka2_micro$genes <- rownames(Scr2vsSka2_micro)
Scr2vsSka2_micro$genes <- as.character(ensembl2gene(Scr2vsSka2_micro$genes))
write.csv(Scr2vsSka2_micro, "/PHShome/je637/RNAseq/RNAseq_Ska2/output/all_genes_2weeks.csv")

Scr2vsSka2_micro <- filter(Scr2vsSka2_micro, padj <= 0.05)
length(rownames(Scr2vsSka2_micro)) 

write.csv(Scr2vsSka2_micro, "/PHShome/je637/RNAseq/RNAseq_Ska2/output/significant_genes_2weeks.csv")
          
Scr4vsSka4_micro <- results(dds_cell_micro, pAdjustMethod = "bonferroni", contrast=c("condition", "SKA2_KD_week4", "SCR_week4"))
Scr4vsSka4_micro <- as.data.frame(Scr4vsSka4_micro)
Scr4vsSka4_micro$genes <- rownames(Scr4vsSka4_micro)
Scr4vsSka4_micro$genes <- as.character(ensembl2gene(Scr4vsSka4_micro$genes))
write.csv(Scr4vsSka4_micro, "/PHShome/je637/RNAseq/RNAseq_Ska2/output/all_genes_4weeks.csv")

Scr4vsSka4_micro <- filter(Scr4vsSka4_micro, padj <= 0.05)
length(rownames(Scr4vsSka4_micro)) 

write.csv(Scr4vsSka4_micro, "/PHShome/je637/RNAseq/RNAseq_Ska2/output/significant_genes_4weeks.csv")


# --------

# Check if normalization went well 

# --------
## PCA after normalization
norm.data <-counts(dds_cell_micro, normalized=TRUE)
pd <- metadata[,c(3:20)] # check the colnames I need
log_data <- log2(norm.data + 1)
custom.PCA(beta = log_data, pd = pd, plot.title = "PCA after normalization")

# --------

# Check for co-linearity
# --------
library(car)
#define multiple linear regression model
model <- lm(num_condition ~ libbatch + prop_ependyma + prop_neuron + prop_oligo + 
              prop_endothelial + prop_neurogenesis + prop_astrocyte + 
              prop_polydendrocyte + prop_mural + prop_fibroblast, data=metadata[,c(3:7,9:19)])
#calculate the VIF for each predictor variable in the model
vif(model)
##VIF equal to 1 = variables are not correlated
##VIF between 1 and 5 = variables are moderately correlated 
##VIF greater than 5 = variables are highly correlated
# There are quite some scores that reach the level of above 5. Therefore there is collinearity within the data

data_x <-  metadata[,c(5:7,9:19)]  # independent variables 
var <- cor(data_x)  # independent variables correlation matrix 
# There we can see that totalng and rnaconc are highly correlated with one another

# --------

# Variance stabilization transformation
# --------
## This is used to reduce type 1 error
vsd <- vst(object=dds_cell_micro,blind=FALSE)
vst_values <- assay(vsd)

custom.PCA(beta = vst_values, pd = coldata, plot.title = "test")

pcaData <- plotPCA.deseq2_modified(vsd, 
                                   ntop = 5000,
                                   intgroup="condition", 
                                   returnData=TRUE)
pca_PCs = pcaData[[1]]
pca_origin = pcaData[[2]]
percentVar <- round(100 * attr(pcaData[[1]], "percentVar"))
res.ind<-get_pca_ind(pca_origin)
fviz_eig(pca_origin) # output is a scree plot for the number of PC's against the percentage of variances that is explained by the PC

# correlation plot
cor_mat <- cbind(pca_PCs[,1:5],
                 metadata[,c(5:7,9:19)])
cor_mat<-sapply(cor_mat, as.double)
rownames(cor_mat)<-rownames(pca_PCs)
all(!is.na(cor_mat))
M <- cor(cor_mat)
res1 <- cor.mtest(cor_mat, conf.level = .95)

corrplot(M, method = "ellipse",type = 'upper',
         tl.pos = 'tp', tl.srt = 45, tl.col = 'black',
         p.mat = res1$p, sig.level = .05, insig = 'blank')
corrplot(M, add = TRUE,method = "number",type = 'lower',
         tl.pos = 'n', cl.pos = 'n',
         p.mat = res1$p, sig.level = .05, insig = 'blank')

pdf(file = paste0(output_path, "SF5.pdf"), width = 15, height = 10)
corrplot(M, method = "ellipse",type = 'upper',
         tl.pos = 'tp', tl.srt = 45, tl.col = 'black',
         p.mat = res1$p, sig.level = .05, insig = 'blank')
corrplot(M, add = TRUE,method = "number",type = 'lower',
         tl.pos = 'n', cl.pos = 'n',
         p.mat = res1$p, sig.level = .05, insig = 'blank')
dev.off()

# -------
# GO and KEGG enrichment analysis
all_genes <- rownames(filtered_data)


## model including cell population of microglia/neurons

sig_week2 <- rownames(Scr2vsSka2_micro)
sig_week4 <- rownames(Scr4vsSka4_micro)
## GO analysis 
enrich_2week <- enrichGO(
  gene = sig_week2,
  OrgDb = "org.Mm.eg.db",
  keyType = "ENSEMBL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "bonferroni",
  universe = all_genes,
  qvalueCutoff = 0.05,
  readable = FALSE)
length(which(enrich_2week@result$p.adjust <= 0.05))
enrich_2week_micro <- enrich_2week@result[c(1:318),]

enrich_4week <- enrichGO(
  gene = sig_week4,
  OrgDb = "org.Mm.eg.db",
  keyType = "ENSEMBL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "bonferroni",
  universe = all_genes,
  qvalueCutoff = 0.05,
  readable = FALSE)
x <- length(which(enrich_4week@result$p.adjust <= 0.05))
enrich_4week_micro <- enrich_4week@result[c(1:x),]

# KEGG analysis
search_kegg_organism('mmu', by='kegg_code')
mmus <- search_kegg_organism("Mus musculus", by="scientific_name")

# Convert Ensembl ID to NCBI ID
ensembl=useMart("ensembl")
ensembl = useDataset("mmusculus_gene_ensembl", mart = ensembl)


# KEGG analysis

# KEGG analysis with microglia in model
week2_KEGG_micro <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
                    filters = "ensembl_gene_id",
                    values = rownames(Scr2vsSka2_micro),
                    mart = ensembl)

week4_KEGG_micro <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
                    filters = "ensembl_gene_id",
                    values = rownames(Scr4vsSka4_micro),
                    mart = ensembl)
all_genes_KEGG <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
                        filters = "ensembl_gene_id",
                        values = all_genes,
                        mart = ensembl)


res2_kegg_micro <- KEGG_enrichment(week2_KEGG_micro$ensembl_gene_id, all_genes_KEGG$ensembl_gene_id)
res4_kegg_micro <- KEGG_enrichment(week4_KEGG_micro$ensembl_gene_id, all_genes_KEGG$ensembl_gene_id)

# Calculate Z-scores

# Zscores calculation for enrichment analysis GO
# ----------
# What I need to do is rerun the enrichment analysis with redable = FALSE

# week 2 with microglia prop
bubble_2_GO_micro <- enrich_2week_micro[,c(1:3,9,7)]
colnames(bubble_2_GO_micro) <- c("Category", "ID", "Term", "Genes", "adj_pval")
bubble_2_GO_micro$Genes <- gsub("/", ",", bubble_2_GO_micro$Genes, fixed = TRUE)
bubble_2_GO_micro$Genes <- toupper(bubble_2_GO_micro$Genes)
bubble_2_GO_micro$Genes <- as.factor(bubble_2_GO_micro$Genes)

bubble_2_sig_micro <- Scr2vsSka2_micro[,c(7,2,1,5,6)]
bubble_2_sig_micro$genes <- rownames(bubble_2_sig_micro)
colnames(bubble_2_sig_micro) <- c("ID", "logFC", "AveExpr", "P.value", "adj.P.Val")

## zscores week 2
zscore_week2_micro <- circle_dat(bubble_2_GO_micro, bubble_2_sig_micro)
zscore_week2_micro <- subset(zscore_week2_micro, !duplicated(zscore_week2_micro$ID))

# week 4
bubble_4_GO_micro <- enrich_4week_micro[,c(1:3,9,7)]
colnames(bubble_4_GO_micro) <- c("Category", "ID", "Term", "Genes", "adj_pval")
bubble_4_GO_micro$Genes <- gsub("/", ",", bubble_4_GO_micro$Genes, fixed = TRUE)
bubble_4_GO_micro$Genes <- toupper(bubble_4_GO_micro$Genes)
bubble_4_GO_micro$Genes <- as.factor(bubble_4_GO_micro$Genes)


bubble_4_sig_micro <- Scr4vsSka4_micro[,c(7,2,1,5,6)]
bubble_4_sig_micro$genes <- rownames(bubble_4_sig_micro)
colnames(bubble_4_sig_micro) <- c("ID", "logFC", "AveExpr", "P.value", "adj.P.Val")

## zscores week 4
zscore_week4_micro <- circle_dat(bubble_4_GO_micro, bubble_4_sig_micro)
zscore_week4_micro <- subset(zscore_week4_micro, !duplicated(zscore_week4_micro$ID))

zscores <- list(zscore_week2_micro, zscore_week4_micro)
names(zscores) <- c("zscore_week2_with_microglia","zscore_week4_with_microglia")
for(i in names(zscores)){
  if(i == "zscore_week2_with_microglia"){
    n = i
    wb <- createWorkbook()
    addWorksheet(wb, n)
    writeData(wb, n, zscores[[i]])
    saveWorkbook(wb, file = paste0(output_path, "zscores_enrichment_analysis.xlsx"), overwrite = TRUE)
  } else {
    n = i
    wb <- loadWorkbook(file = paste0(output_path, "zscores_enrichment_analysis.xlsx"))
    addWorksheet(wb, n)
    writeData(wb, n, zscores[[i]], startRow = 1, startCol = 1, colNames = TRUE)
    saveWorkbook(wb, file = paste0(output_path, "zscores_enrichment_analysis.xlsx"), overwrite = TRUE)
    
  }
}
# ----------

# Zscores calculation for enrichment analysis KEGG
# -----

# week 2  microglia proportion
bubble_2_KEGG_micro <- res2_kegg_micro[,c(1,2,6,8)]

colnames(bubble_2_KEGG_micro) <- c("Category", "Term", "adj_pval","Genes")
bubble_2_KEGG_micro$Genes <- gsub("/", ",", bubble_2_KEGG_micro$Genes, fixed = TRUE)
bubble_2_KEGG_micro$Genes <- toupper(bubble_2_KEGG_micro$Genes)
bubble_2_KEGG_micro$Genes <- as.factor(bubble_2_KEGG_micro$Genes)


bubble_2_sig_micro <- Scr2vsSka2_micro[,c(7,2,1,5,6)]
bubble_2_sig_micro$genes <- rownames(bubble_2_sig_micro)
colnames(bubble_2_sig_micro) <- c("ID", "logFC", "AveExpr", "P.value", "adj.P.Val")
zscore_week2_KEGG_micro <- GOplot::circle_dat(bubble_2_KEGG_micro, bubble_2_sig_micro)

zscore_week2_KEGG_micro <- subset(zscore_week2_KEGG_micro, !duplicated(zscore_week2_KEGG_micro$term))


# -----

save.image("/PHShome/je637/RNAseq/RNAseq_Ska2/workenvironment/DEG_analysis_cellprop.RDS")


