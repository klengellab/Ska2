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
options(stringsAsFactors = FALSE)
set.seed(1234)

# Functions
source("/PHShome/je637/gitlab/general_functions/plot/customPCA.R")
ensembl2gene <- function(ensembl_name){
  if(length(ensembl_name > 1)){
    gene <- list()
    for(i in ensembl_name){
      ens_id <- genes %>% dplyr::filter(ensembl_gene_id == i) %>%
        dplyr::select(external_gene_name) %>% pull()
      if(length(ens_id) == 0){
        gene <- append(gene, i)
      } else if(ens_id == ""){
        gene <- append(gene, i)
      } else if(length(ens_id != 0)){
        gene <- append(gene, ens_id)
      } else {
        gene <- append(gene, i)
      }
    }
  } else {
    ens_id <- genes %>% dplyr::filter(ensembl_gene_id == ensembl_name) %>%
      dplyr::select(external_gene_name) %>% pull()
    if(length(ens_id) == 0){
      gene <- append(gene, i)
    } else if(ens_id == ""){
      gene <- append(gene, i)
    } else if(length(ens_id != 0)){
      gene <- append(gene, ens_id)
    } else {
      gene <- append(gene, i)
    }
  }
  return(gene)
}
genes <- read.csv("/PHShome/je637/general/tables/ensembl_w_description.mouse.csv")
biomart_ensembl_to_gene <- function(gene){
  ensembl=useMart("ensembl")
  ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
  names <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                 filters = 'ensembl_gene_id', 
                 values = gene, 
                 mart = ensembl)
  return(names)
}
Volcanoplot <- function(data, title){
  # This function will plot a volcano plot of the DEGs with nominal and FDR significant
  # genes
  # Data must be a dataframe
  # title is character string of the title of the volcanoplot
  
  data$dif <- "NS"
  data$dif[data[,2] > 0 & data$pvalue < 0.05] <- "UP-nom"
  data$dif[data[,2] < 0 & data$pvalue < 0.05] <- "DOWN-nom"
  data$dif[data[,2] > 0 & data$padj < 0.05] <- "UP"
  data$dif[data[,2] < 0 & data$padj < 0.05] <- "DOWN"
  mycolors <- c("#4169E1","#7d9bf5", "#AFADB3", "#c97171","#B22222")
  names(mycolors) <- c("DOWN","DOWN-nom", "NS", "UP-nom", "UP")
  data$symbol <- rownames(data)
  data$delabel <- NA
  data$delabel[data$dif %in% c("UP", "DOWN")] <- data$symbol[data$dif %in% c("UP", "DOWN")]
  test <- biomart_ensembl_to_gene(data$delabel)
  x <- which(data$delabel %in% test$ensembl_gene_id)
  data[x,10] <- test$external_gene_name
  ggrepel.max.overlaps = Inf
  p <- ggplot(data = data, aes(x=data[,2], y=-log10(pvalue), col=dif)) + geom_point() + theme_classic()
  p2 <- p + scale_color_manual(values = mycolors) 
  p3 <- p2 + labs(title = title) + xlab("Log2FC") + ylab("-log10 p.value")
  return(p3)
}

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
pd <- metadata[,c(3:8,9:19)]
log_data <- log2(filtered_data + 1)
custom.PCA(beta = log_data, pd = pd, plot.title = "PCA before normalization")



# --------

# Differential expression analysis with Deseq2

# --------
# Differential expression analysis with Deseq2 including the cell type proportions
coldata <- metadata[,c(1,3:8,9:19)]
coldata$libbatch <- as.factor(coldata$libbatch)
rownames(coldata) <- coldata$samplename

## excluding cell proportions
#-------
all(rownames(coldata) == colnames(filtered_data)) # Make sure the output is true
dds <- DESeqDataSetFromMatrix(countData = filtered_data, colData = coldata, 
                              design = ~ libbatch + condition)
dds$condition <- factor(dds$condition, levels = c("SCR_week2", "SKA2_KD_week2", "SCR_week4", "SKA2_KD_week4"))
dds <- DESeq(dds)
resultsNames(dds)
dds0 <- dds # for later use with surrogate variable analysis
#-------

## including all cell proportions
#-------
dds_cellprop <- DESeqDataSetFromMatrix(countData = filtered_data, colData = coldata,
                                       design = ~ libbatch + prop_ependyma + prop_oligo + 
                                         prop_endothelial + prop_neurogenesis + prop_astrocyte + 
                                         prop_polydendrocyte + prop_mural + prop_fibroblast + prop_microglia + condition)
dds_cellprop$condition <- factor(dds_cellprop$condition, levels = c("SCR_week2", "SKA2_KD_week2", "SCR_week4", "SKA2_KD_week4"))
dds_cellprop <- DESeq(dds_cellprop)
resultsNames(dds_cellprop)
#-------

## including microglia population
#-------
dds_cell_micro <- DESeqDataSetFromMatrix(countData = filtered_data, colData = coldata,
                                       design = ~ libbatch + prop_microglia + condition)
dds_cell_micro$condition <- factor(dds_cell_micro$condition, levels = c("SCR_week2", "SKA2_KD_week2", "SCR_week4", "SKA2_KD_week4"))
dds_cell$condition <- factor(dds_cell$condition, levels = c("SCR_week2", "SKA2_KD_week2", "SCR_week4", "SKA2_KD_week4"))
dds_cell_micro <- DESeq(dds_cell_micro)
resultsNames(dds_cell_micro)
#-------

## including microglia and mural cell population
#--------
dds_cell_micro_mural <- DESeqDataSetFromMatrix(countData = filtered_data, colData = coldata,
                                         design = ~ libbatch + prop_microglia + prop_mural + condition)

dds_cell_micro_mural$condition <- factor(dds_cell_micro_mural$condition, levels = c("SCR_week2", "SKA2_KD_week2", "SCR_week4", "SKA2_KD_week4"))
dds_cell_micro_mural <- DESeq(dds_cell_micro_mural)
resultsNames(dds_cell_micro_mural)
                           
# Check if normalization went well 

# --------
## PCA after normalization
norm.data <-counts(dds_cell_micro, normalized=TRUE)
pd <- metadata[,c(3:8,9:19)] # check the colnames I need
log_data <- log2(norm.data + 1)
custom.PCA(beta = log_data, pd = pd, plot.title = "PCA after normalization")



## Results 2 weeks

# without cell type proportions
Scr2vsSka2 <- results(dds, pAdjustMethod = "bonferroni", contrast=c("condition", "SKA2_KD_week2", "SCR_week2"))
sigScr2vsSka2 <- as.data.frame(Scr2vsSka2)
sigScr2vsSka2$genes <- rownames(Scr2vsSka2)
sig2 <- filter(sigScr2vsSka2, padj <= 0.05)
length(rownames(sig2)) #3482

# with all cell proporotions
Scr2vsSka2_prop <- results(dds_cellprop, pAdjustMethod = "bonferroni", contrast=c("condition", "SKA2_KD_week2", "SCR_week2"))
sigScr2vsSka2_prop <- as.data.frame(Scr2vsSka2_prop)
sigScr2vsSka2_prop$genes <- rownames(Scr2vsSka2_prop)
sig2_prop <- filter(sigScr2vsSka2_prop, padj <= 0.05)
length(rownames(sig2_prop)) #256

# only with microglia proportions included in the model
Scr2vsSka2_micro <- results(dds_cell_micro, pAdjustMethod = "bonferroni", contrast=c("condition", "SKA2_KD_week2", "SCR_week2"))
Scr2vsSka2_micro <- as.data.frame(Scr2vsSka2_micro)
Scr2vsSka2_micro$genes <- rownames(Scr2vsSka2_micro)
write.csv(Scr2vsSka2_micro, "/PHShome/je637/RNAseq/RNAseq_Ska2/output/all_genes_2weeks_V2.csv")
Scr2vsSka2_micro <- filter(Scr2vsSka2_micro, padj <= 0.05)
length(rownames(Scr2vsSka2_micro)) #1367
write.csv(Scr2vsSka2_micro, "/PHShome/je637/RNAseq/RNAseq_Ska2/output/significant_genes_2weeks_V2.csv")


# including mural and microglia in the model
Scr2vsSka2_mm <- results(dds_cell_micro_mural, pAdjustMethod = "bonferroni", contrast=c("condition", "SKA2_KD_week2", "SCR_week2"))
Scr2vsSka2_mm <- as.data.frame(Scr2vsSka2_mm)
Scr2vsSka2_mm$genes <- rownames(Scr2vsSka2_mm)
Scr2vsSka2_mm_sig <- filter(Scr2vsSka2_mm, padj <= 0.05)
length(rownames(Scr2vsSka2_mm_sig)) #1012


# Venn diagram of the overlap between 2 week significant genes with and without adjusting for cell proportions
venn.diagram(x = list(rownames(sig2), rownames(Scr2vsSka2_micro), rownames(sig2_prop), rownames(Scr2vsSka2_mm)),
             category.names = c("DEGs week 2","DEGs with micro", "DEGs with all celltypes", "DEGs with mm"),
             filename = "week2_DEG_overlap.png")

#Volcanoplots

# without cell type composition
library(ggrepel)

pdf(paste0(output_path, "volcanoplots_without_cell_type.pdf"))
t <- Volcanoplot(as.data.frame(Scr2vsSka2), "control vs. Ska2 at 2 weeks")
plot(t)
dev.off()

pdf(paste0(output_path, "volcanoplots_with_microglia.pdf"))
t <- Volcanoplot(as.data.frame(Scr2vsSka2_micro), "control vs. Ska2 at 2 weeks with microglia prop")
plot(t)
dev.off()

pdf(paste0(output_path, "volcanoplots_with_mm.pdf"))
t <- Volcanoplot(as.data.frame(Scr2vsSka2_mm), "control vs. Ska2 at 2 weeks with micro and mural")
plot(t)
dev.off()

pdf(paste0(output_path, "volcanoplots_with_all.pdf"))
t <- Volcanoplot(as.data.frame(Scr2vsSka2_prop), "control vs. Ska2 at 2 weeks with all cells")
plot(t)
dev.off()

# Results differential expression analysis week 4

# without celltypes in the model
Scr4vsSka4 <- results(dds, pAdjustMethod = "bonferroni", contrast=c("condition", "SKA2_KD_week4", "SCR_week4"))
sigScr4vsSka4 <- as.data.frame(Scr4vsSka4)
sigScr4vsSka4$genes <- rownames(Scr4vsSka4)
sigScr4 <- filter(sigScr4vsSka4, padj <= 0.05) #5573

# with all celltypes in the model
Scr4vsSka4_prop <- results(dds_cellprop, pAdjustMethod = "bonferroni", contrast=c("condition", "SKA2_KD_week4", "SCR_week4"))
sigScr4vsSka4_prop <- as.data.frame(Scr4vsSka4_prop)
sigScr4vsSka4_prop$genes <- rownames(Scr4vsSka4_prop)
sigScr4_prop <- filter(sigScr4vsSka4_prop, padj <= 0.05) #78

# with microglia cells in the model
Scr4vsSka4_micro <- results(dds_cell_micro, pAdjustMethod = "bonferroni", contrast=c("condition", "SKA2_KD_week4", "SCR_week4"))
Scr4vsSka4_micro <- as.data.frame(Scr4vsSka4_micro)
Scr4vsSka4_micro$genes <- rownames(Scr4vsSka4_micro)
Scr4vsSka4_micro$genes <- as.character(ensembl2gene(Scr4vsSka4_micro$genes))
write.csv(Scr4vsSka4_micro, "/PHShome/je637/RNAseq/RNAseq_Ska2/output/all_genes_4weeks_V2.csv")
Scr4vsSka4_micro <- filter(Scr4vsSka4_micro, padj <= 0.05)
length(rownames(Scr4vsSka4_micro)) #601
write.csv(Scr4vsSka4_micro, "/PHShome/je637/RNAseq/RNAseq_Ska2/output/significant_genes_4weeks_V2.csv")

# model with microglia and mural cells
Scr4vsSka4_mm <- results(dds_cell_micro_mural, pAdjustMethod = "bonferroni", contrast=c("condition", "SKA2_KD_week4", "SCR_week4"))
Scr4vsSka4_mm <- as.data.frame(Scr4vsSka4_mm)
Scr4vsSka4_mm$genes <- rownames(Scr4vsSka4_mm)
Scr4vsSka4_mm_sig <- filter(Scr4vsSka4_mm, padj <= 0.05)
length(rownames(Scr4vsSka4_mm_sig)) #646

# Venn diagram of the overlap between 2 week significant genes with and without adjusting for cell proportions
venn.diagram(x = list(rownames(sigScr4), rownames(Scr4vsSka4_micro), rownames(sigScr4_prop), rownames(Scr4vsSka4_mm_sig)),
             category.names = c("DEGs week 4", "DEGs week 4 with microglia", "DEGs week 4 with all celltypes", "DEGs week 4 with mm"),
             filename = "week4_DEG_overlap.png")

# Volcanoplots 4 weeks
pdf(paste0(output_path, "volcanoplots_without_cell_type_week4.pdf"))
t <- Volcanoplot(as.data.frame(Scr4vsSka4), "control vs. Ska2 at 4 weeks")
plot(t)
dev.off()

pdf(paste0(output_path, "volcanoplots_with_microglia_week4.pdf"))
t <- Volcanoplot(as.data.frame(Scr4vsSka4_micro), "control vs. Ska2 at 4 weeks with microglia prop")
plot(t)
dev.off()

pdf(paste0(output_path, "volcanoplots_with_mm_week4.pdf"))
t <- Volcanoplot(as.data.frame(Scr4vsSka4_mm), "control vs. Ska2 at 4 weeks with micro and mural")
plot(t)
dev.off()

pdf(paste0(output_path, "volcanoplots_with_all_week4.pdf"))
t <- Volcanoplot(as.data.frame(Scr4vsSka4_prop), "control vs. Ska2 at 4 weeks with all celtypes")
plot(t)
dev.off()

# Check if normalization went well 
# --------
## PCA after normalization
norm.data <-counts(dds, normalized=TRUE)
norm.data_cell <- counts(dds_cellprop, normalized = TRUE)
barplot(log2(norm.data + 1))
all(norm.data == norm.data_cell) # This should be true otherwise something is going wrong in your Deseq2 normalization
# since adding covariates to the data shouldn't matter on the normalization of the data.
rm(norm.data_cell)

pd <- metadata[,c(3:8,9:19)] # check the colnames I need
log_data <- log2(norm.data + 1)
custom.PCA(beta = log_data, pd = pd, plot.title = "PCA after normalization")
# --------

# Check for co-linearity
# --------
library(car)
#define multiple linear regression model
model <- lm(num_condition ~ libbatch + prop_ependyma + prop_microglia + prop_oligo + 
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
source("/PHShome/je637/gitlab/rna-seq/Jonathan Human AD/plot_PCA_deseq2_modified.R")
vsd <- vst(object=dds,blind=FALSE)
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

pdf(file = "/PHShome/je637/RNAseq/RNAseq_Ska2/output/correlation_plot_V2.pdf", width = 15, height = 10)
corrplot(M, method = "ellipse",type = 'upper',
         tl.pos = 'tp', tl.srt = 45, tl.col = 'black',
         p.mat = res1$p, sig.level = .05, insig = 'blank')
corrplot(M, add = TRUE,method = "number",type = 'lower',
         tl.pos = 'n', cl.pos = 'n',
         p.mat = res1$p, sig.level = .05, insig = 'blank')
dev.off()
# Here we see that microglia and neuronal cell populations are highly correlated with one another. PC4 captures the library batch
# --------


# Enrichment analysis on the model including microglia
# --------
# GO and KEGG enrichment analysis
all_genes <- rownames(filtered_data)


## model including cell population of microglia/neurons
sig_week2 <- rownames(sig2)
sig_week4 <- rownames(sigScr4)
sig_week2_micro <- rownames(Scr2vsSka2_micro)
sig_week4_micro <- rownames(Scr4vsSka4_micro)
sig_week4_mm <- rownames(Scr4vsSka4_mm_sig)
sig_week2_mm <- rownames(Scr2vsSka2_mm_sig)
sig_week2_all <- rownames(sig2_prop)
sig_week4_all <- rownames(sigScr4_prop)

## GO analysis with all cell types
enrich_2week_all <- enrichGO(
  gene = sig_week2_all,
  OrgDb = "org.Mm.eg.db",
  keyType = "ENSEMBL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "bonferroni",
  universe = all_genes,
  qvalueCutoff = 0.05,
  readable = FALSE)
x <- length(which(enrich_2week_all@result$p.adjust <= 0.05))
enrich_2week_all <- enrich_2week_all@result[1:x,]
dim(enrich_2week_all) #166 10

## GO analysis with microglia and mural
enrich_2week_mm <- enrichGO(
  gene = sig_week2_mm,
  OrgDb = "org.Mm.eg.db",
  keyType = "ENSEMBL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "bonferroni",
  universe = all_genes,
  qvalueCutoff = 0.05,
  readable = FALSE)
x <- length(which(enrich_2week_mm@result$p.adjust <= 0.05))
enrich_2week_mm <- enrich_2week_mm@result[1:x,]
dim(enrich_2week_mm) #240 10

## GO analysis without cell type
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
x <- length(which(enrich_2week@result$p.adjust <= 0.05))
enrich_2week <- enrich_2week@result[1:x,]
dim(enrich_2week) #342 10

## GO with microglia cells
enrich_2week_micro <- enrichGO(
  gene = sig_week2_micro,
  OrgDb = "org.Mm.eg.db",
  keyType = "ENSEMBL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "bonferroni",
  universe = all_genes,
  qvalueCutoff = 0.05,
  readable = FALSE)
x <- length(which(enrich_2week_micro@result$p.adjust <= 0.05))
enrich_2week_micro <- enrich_2week_micro@result[1:x,]
dim(enrich_2week_micro) #254 10

# week 4
sig_week4_mm <- rownames(Scr4vsSka4_mm_sig)
sig_week2_mm <- rownames(Scr2vsSka2_mm_sig)
sig_week2_all <- rownames(sig2_prop)
sig_week4_all <- rownames(sigScr4_prop)

# with microglia and mural cells
enrich_4week_mm <- enrichGO(
  gene = sig_week4_mm,
  OrgDb = "org.Mm.eg.db",
  keyType = "ENSEMBL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "bonferroni",
  universe = all_genes,
  qvalueCutoff = 0.05,
  readable = FALSE)
x <- length(which(enrich_4week_mm@result$p.adjust <= 0.05))
enrich_4week_mm <- enrich_4week_mm@result[c(1:x),]
dim(enrich_4week_mm) # 9 10

# with all cell types 
enrich_4week_all <- enrichGO(
  gene = sig_week4_all,
  OrgDb = "org.Mm.eg.db",
  keyType = "ENSEMBL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "bonferroni",
  universe = all_genes,
  qvalueCutoff = 0.05,
  readable = FALSE)
x <- length(which(enrich_4week_all@result$p.adjust <= 0.05))
enrich_4week_all <- enrich_4week_all@result[c(1:x),]
dim(enrich_4week_all) # 1 10

# without cell types
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
enrich_4week <- enrich_4week@result[c(1:x),]
dim(enrich_4week) # 153 10

# with microglia cell type
enrich_4week_micro <- enrichGO(
  gene = sig_week4_micro,
  OrgDb = "org.Mm.eg.db",
  keyType = "ENSEMBL",
  ont = "ALL",
  pvalueCutoff = 0.05,
  pAdjustMethod = "bonferroni",
  universe = all_genes,
  qvalueCutoff = 0.05,
  readable = FALSE)
x <- length(which(enrich_4week_micro@result$p.adjust <= 0.05))
enrich_4week_micro <- enrich_4week_micro@result[1:x,]
dim(enrich_4week_micro) # 3 10

# KEGG analysis
search_kegg_organism('mmu', by='kegg_code')
mmus <- search_kegg_organism("Mus musculus", by="scientific_name")

## model without cell type proportions

# Convert Ensembl ID to NCBI ID
ensembl=useMart("ensembl")
ensembl = useDataset("mmusculus_gene_ensembl", mart = ensembl)

# I can't use the same KEGG analysis as previously. However that still relied on the clusterprofiler
# Therefore I'm going to use msigdbr

# KEGG analysis

KEGG_enrichment <- function(genes, universe){
  # Enrichment analysis with msigdbr/clusterProfiler
  # genes is a character vector of external gene names
  # univere is a character vector of external gene names
  # Databases are from GSEA
  message("KEGG analysis")
  KEGG_datasets <- msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP:KEGG")
  enriched_KEGG <- as.data.frame(enricher(genes, universe = universe, pvalueCutoff = 0.05, 
                                          qvalueCutoff = 0.05, pAdjustMethod = "BH", TERM2GENE =
                                            KEGG_datasets[,c("gs_name", "ensembl_gene")]))
  enriched_KEGG$Description <- str_sub(enriched_KEGG$ID,6) 
  enriched_KEGG$Description <- str_replace_all(enriched_KEGG$Description, "_", " ")
  
  return(enriched_KEGG)
}

# KEGG analysis without cell type in model
week2_KEGG <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
                    filters = "ensembl_gene_id",
                    values = rownames(sig2),
                    mart = ensembl)
week2_KEGG_micro <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
                          filters = "ensembl_gene_id",
                          values = rownames(Scr2vsSka2_micro),
                          mart = ensembl)

week2_KEGG_mm <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
                       filters = "ensembl_gene_id",
                       values = rownames(Scr2vsSka2_mm_sig),
                       mart = ensembl)

week2_KEGG_all <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
                        filters = "ensembl_gene_id",
                        values = rownames(sig2_prop),
                        mart = ensembl)

week4_KEGG <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
                    filters = "ensembl_gene_id",
                    values = rownames(sigScr4),
                    mart = ensembl)
week4_KEGG_micro <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
                          filters = "ensembl_gene_id",
                          values = rownames(Scr4vsSka4_micro),
                          mart = ensembl)

week4_KEGG_mm <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
                       filters = "ensembl_gene_id",
                       values = rownames(Scr4vsSka4_mm_sig),
                       mart = ensembl)

week4_KEGG_all <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
                        filters = "ensembl_gene_id",
                        values = rownames(sigScr4_prop),
                        mart = ensembl)

all_genes_KEGG <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
                        filters = "ensembl_gene_id",
                        values = all_genes,
                        mart = ensembl)

res2_kegg <- KEGG_enrichment(week2_KEGG$ensembl_gene_id, all_genes_KEGG$ensembl_gene_id)
res2_kegg_micro <- KEGG_enrichment(week2_KEGG_micro$ensembl_gene_id, all_genes_KEGG$ensembl_gene_id)
res2_kegg_all <- KEGG_enrichment(week2_KEGG_all$ensembl_gene_id, all_genes_KEGG$ensembl_gene_id)
res2_kegg_mm <- KEGG_enrichment(week2_KEGG_mm$ensembl_gene_id, all_genes_KEGG$ensembl_gene_id)



res4_kegg <- KEGG_enrichment(week4_KEGG$ensembl_gene_id, all_genes_KEGG$ensembl_gene_id)
res4_kegg_micro <- KEGG_enrichment(week4_KEGG_micro$ensembl_gene_id, all_genes_KEGG$ensembl_gene_id)
res4_kegg_all <- KEGG_enrichment(week4_KEGG_all$ensembl_gene_id, all_genes_KEGG$ensembl_gene_id)
res4_kegg_mm <- KEGG_enrichment(week4_KEGG_mm$ensembl_gene_id, all_genes_KEGG$ensembl_gene_id)



# Calculate Z-scores
# This is not yet working nicely, check the preprocessing_Final.R file on your own home directory
# as reference to what you did earlier

# output to excel sheets

week2 <- list("enrichment GO" = enrich_2week, "enrichment GO with microglia" = enrich_2week_micro,
              "GO enrichment with mm" = enrich_2week_mm ,
              " enrichment GO with all" = enrich_2week_all, "enrichment KEGG" = res2_kegg,
              "enrichment KEGG with microglia" = res2_kegg_micro,
              "enrichment KEGG with mm" = res2_kegg_mm, "enrichment KEGG with all" = res2_kegg_all)
week4 <- list("enrichment GO" = enrich_4week, "enrichment GO with microglia" = enrich_4week_micro,
              "GO enrichment with mm" = enrich_4week_mm ,
              " enrichment GO with all" = enrich_4week_all, "enrichment KEGG" = res4_kegg,
              "enrichment KEGG with microglia" = res4_kegg_micro,
              "enrichment KEGG with mm" = res4_kegg_mm, "enrichment KEGG with all" = res4_kegg_all)

# Add the KEGG analysis to week 2 and 4 as well
for(i in names(week2)){
  if(i == "enrichment GO"){
    n = i
    wb <- createWorkbook()
    addWorksheet(wb, n)
    writeData(wb, n, week2[[i]])
    saveWorkbook(wb, file = "/PHShome/je637/RNAseq/RNAseq_Ska2/output/enrichment_analysis_week2_V2.xlsx", overwrite = TRUE)
  } else {
    n = i
    wb <- loadWorkbook(file = "/PHShome/je637/RNAseq/RNAseq_Ska2/output/enrichment_analysis_week2_V2.xlsx")
    addWorksheet(wb, n)
    writeData(wb, n, week2[[i]], startRow = 1, startCol = 1, colNames = TRUE)
    saveWorkbook(wb, file = "/PHShome/je637/RNAseq/RNAseq_Ska2/output/enrichment_analysis_week2_V2.xlsx", overwrite = TRUE)

  }
}

for(i in names(week4)){
  if(i == "enrichment GO"){
    n = i
    wb <- createWorkbook()
    addWorksheet(wb, n)
    writeData(wb, n, week4[[i]])
    saveWorkbook(wb, file = "/PHShome/je637/RNAseq/RNAseq_Ska2/output/enrichment_analysis_week4_V2.xlsx", overwrite = TRUE)
  } else {
    n = i
    wb <- loadWorkbook(file = "/PHShome/je637/RNAseq/RNAseq_Ska2/output/enrichment_analysis_week4_V2.xlsx")
    addWorksheet(wb, n)
    writeData(wb, n, week4[[i]], startRow = 1, startCol = 1, colNames = TRUE)
    saveWorkbook(wb, file = "/PHShome/je637/RNAseq/RNAseq_Ska2/output/enrichment_analysis_week4_V2.xlsx", overwrite = TRUE)
    
  }
}

# ------

# Correlation analysis
# -------- 
# condition vs. microglia prop/neuron prop
cor(metadata$num_condition, metadata$prop_microglia) #0.1169199
cor(metadata$num_condition, metadata$prop_neuron) #-0.1428665

# microglia vs. neuron
cor(metadata$prop_microglia, metadata$prop_neuron) #-0.9726371
# -------- 




