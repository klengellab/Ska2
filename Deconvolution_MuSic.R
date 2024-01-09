# Deconvolution of Ska2 data with MuSiC

# libraries
library(Biobase)
library(limma)
library(marray)
library(DropSeq.util)
library(biomaRt)
library(MuSiC)
library(ggplot2)
library(nnls)
library(xbioc)
library(reshape)
library(cowplot)
library(DropSeq.util)
library(SingleCellExperiment)
library(reshape2)

# Load in data regarding snRNA-seq data
scRNA_anno = readRDS("/PHShome/je637/Datasets/Dropviz/annotation.BrainCellAtlas_Saunders_version_2018.04.01.RDS")
scRNA_expr = readRDS("/PHShome/je637/Datasets/Dropviz/metacells.BrainCellAtlas_Saunders_version_2018.04.01.RDS")
outcomes_hippo = readRDS("/PHShome/je637/Datasets/Dropviz/F_GRCm38.81.P60Hippocampus.cell_cluster_outcomes.RDS")
dge.path <- "/PHShome/je637/Datasets/Dropviz/F_GRCm38.81.P60Hippocampus.raw.dge.txt.gz"

# Open data
dge <- loadSparseDge(dge.path)
Expr_hippo <- as.matrix(dge)
rm(dge)

# Select for only samples that belong to the Hippocampus region annotation file (phenodata)
HC <- which(scRNA_anno$tissue == "HC")
scRNA_hippo_anno <- scRNA_anno[HC,]
scRNA_hippo_anno <- as.data.frame(scRNA_hippo_anno)

# Select for only samples that belong to the Hippocampus region expression file (expressiondata)
scRNA_hippo_expr <- scRNA_expr[,grepl("HC", colnames(scRNA_expr))]

# Filter out cells that do not contribute e.a. doublets, outliers
out <- which(outcomes_hippo$reason %in% c("curation", "doublet", "min_genes", "not_enough_ics", "outlier", "singleton_cluster", "singleton_subcluster", "small_cell"))
# We are deleting out 81584 cells
outcomes_hippo <- outcomes_hippo[-c(out), ]
outcomes_hippo$samples <- rownames(outcomes_hippo)

# filter out all samples with doublets, outliers in expression data set
all(rownames(outcomes_hippo) %in% colnames(Expr_hippo)) #TRUE
all(rownames(outcomes_hippo) == colnames(Expr_hippo)) #FALSE
reorder <- match(rownames(outcomes_hippo), colnames(Expr_hippo))
Expr_hippo <- Expr_hippo[, reorder]
all(rownames(outcomes_hippo) == colnames(Expr_hippo)) #TRUE
dim(Expr_hippo) #27953 52846

rownames_Expr_hippo <- as.data.frame(rownames(Expr_hippo))
rownames_Expr_hippo$genes <- as.character(rownames_Expr_hippo$`rownames(Expr_hippo)`)

#Convert geneID rownames to Ensembl IDs
ensembl <- useMart("ensembl")
ensembl = useDataset("mmusculus_gene_ensembl", mart = ensembl)
IDs = getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
            filters = 'external_gene_name',
            values = rownames_Expr_hippo$genes,
            mart = ensembl)
dim(IDs) # 24082 2

# The number of genes with external gene names by biomart are less compared to the ones in
# in the single cell dataset available. Therefore we need to delete out the genes that are not found by biomart. Since we need similar gene names for the deconvolution to work.
x <- which(IDs$external_gene_name %in% rownames(Expr_hippo))
IDs <- IDs[x,]
all(IDs$external_gene_name %in% rownames(Expr_hippo)) #TRUE
all(rownames(Expr_hippo) == IDs$external_gene_name) # False
reorder <- match(IDs$external_gene_name, rownames(Expr_hippo))
Expr_hippo <- Expr_hippo[reorder, ] #24081 52846
all(rownames(Expr_hippo) == IDs$external_gene_name)
rownames(Expr_hippo) <- IDs$ensembl_gene_id
message("finished part 7")
outcomes_hippo$cluster <- as.character(outcomes_hippo$cluster)
outcomes_hippo$cluster[outcomes_hippo$cluster %in% as.character(c(1:6, 14))] <- "NEURON"
outcomes_hippo$cluster[outcomes_hippo$cluster == "7"] <- "ASTROCYTE"
outcomes_hippo$cluster[outcomes_hippo$cluster == "10"] <- "MICROGLIA"
outcomes_hippo$cluster[outcomes_hippo$cluster == "8"] <- "OLIGODENDROCYTE"
outcomes_hippo$cluster[outcomes_hippo$cluster == "9"] <- "POLYDENDROCYTE"
outcomes_hippo$cluster[outcomes_hippo$cluster == "15"] <- "ENDOTHELIAL"
outcomes_hippo$cluster[outcomes_hippo$cluster == "17"] <- "FIBROBLAST"
outcomes_hippo$cluster[outcomes_hippo$cluster == "16"] <- "MURAL"
outcomes_hippo$cluster[outcomes_hippo$cluster == "13"] <- "NEUROGENESIS"
outcomes_hippo$cluster[outcomes_hippo$cluster %in% as.character(c(12, 19))] <- "CP"
outcomes_hippo$cluster[outcomes_hippo$cluster == "11"] <- "EPENDYMA"
outcomes_hippo <- outcomes_hippo[,c(1,2,4)] # deleting the column with reason (usually NA but the QC details were stored over here such as doublet etc

# Load in the bulk RNA-seq dataset
data <- read.csv("/PHShome/je637/ska2/data/counts.csv", header = TRUE)
metadata <- read.csv("/PHShome/je637/ska2/data/jakob_ska2.csv", header = TRUE, sep = ",")
rownames(data) <- data$gene
data <- data[,2:33]

# Matching count data and metadata
metadata$samplename <- paste0("X", metadata$samplename)
rownames(metadata) <- metadata$samplename
all(rownames(metadata) %in% colnames(data))
all(rownames(metadata) == colnames(data))
data <- data[, match(rownames(metadata), colnames(data))]
all(rownames(metadata) == colnames(data)) # Outcome should be true
meta <- AnnotatedDataFrame(metadata)
bulk.eset <- ExpressionSet(as.matrix(data), phenoData = meta)
bulk.eset <- as.matrix(bulk.eset)

# Running MuSiC
sce <- SingleCellExperiment(assays = list(counts=Expr_hippo), colData=outcomes_hippo)
Est.prop.Jakob = music_prop(bulk.mtx = bulk.eset, sc.sce = sce, clusters = 'cluster',
                            samples = 'samples', verbose = F)

# Plot the results of MuSiC
Results <- as.data.frame(Est.prop.Jakob[["Est.prop.weighted"]])

m.prop = rbind(melt(Est.prop.Jakob$Est.prop.weighted))
groups <- paste0(metadata$group, "_", metadata$timepoint)
groups <- rep(groups, 10)
m.prop$exp.condition <- groups


Selection<- m.prop[c(33:96, 161:224, 257:320),]
p <- ggplot(Selection, aes(exp.condition, value)) + geom_jitter(aes(colour = exp.condition))
p1 <- p + facet_grid(cols = vars(Var2)) + theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
                                                axis.ticks.x = element_blank()) + labs(y = "Est. prop cells")

pdf("/PHShome/je637/ska2/output/Deconvolution_V2.pdf")
plot(p1)
dev.off()

write.csv(Selection, "/PHShome/je637/ska2/output/proportions_celltypes.csv")
save.image("/PHShome/je637/RNAseq/RNAseq_Ska2/workenvironment/Music_Ska2.RData")
