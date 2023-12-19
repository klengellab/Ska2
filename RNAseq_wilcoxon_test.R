# Wilcoxon Rank sum test to determine the overlap between RNA-seq and secretrome data of the Ska2 study
# Date: 12/09/20
# By: Joy Otten

# reference: https://stat.ethz.ch/R-manual/R-devel/library/stats/html/wilcox.test.html

# libraries:
library(biomaRt)
library(ggplot2)
set.seed(1234)

# Functions
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

# Read in Data
secretome <- read.csv(paste0(output_path, "secretome_JH_2020_dec.csv"),
                      header = TRUE, sep = ",")
genes_2weeks <- read.csv(paste0(output_path, "genes_2weeks.csv"))
genes_4weeks <- read.csv(paste0(output_path, "genes_4weeks.csv"))

# Select only for the significant genes in the secretome data 
secretome_sig <- secretome[,c(7,2)]
x <- which(secretome_sig$P.value < 0.05)
secretome_sig <- secretome_sig[x,]
# In total there are 446 genes significant different (p.value < 0.05) in the secretome data

# Biomart to convert ensemble id to external gene names
gene_names_2weeks <- biomart_ensembl_to_gene(genes_2weeks$genes)
gene_names_4weeks <- biomart_ensembl_to_gene(genes_4weeks$genes)

x <- match(gene_names_2weeks$ensembl_gene_id, genes_2weeks$genes)
y <- match(gene_names_4weeks$ensembl_gene_id, genes_4weeks$genes)

genes_2weeks <- genes_2weeks[x,]
genes_4weeks <- genes_4weeks[y,]

genes_2weeks$gene_names <- gene_names_2weeks$external_gene_name
genes_4weeks$gene_names <- gene_names_4weeks$external_gene_name

# wilcoxon rank sum test
# -----------
length(which(secretome_sig$Gene.names %in% genes_2weeks$gene_names)) #429
length(which(secretome_sig$Gene.names %in% genes_4weeks$gene_names)) #429

length(genes_2weeks$padj) == length(unique(genes_2weeks$padj))
length(unique(genes_2weeks$padj)) #4507
# There are quit some genes with the same padj therefore I am using the nominal p.values
length(unique(genes_2weeks$pvalue)) #20443

length(genes_4weeks$padj) == length(unique(genes_4weeks$padj))
length(unique(genes_4weeks$padj)) #6791
# There are quit some genes with the same padj therefore I am using the nominal p.values
length(unique(genes_2weeks$pvalue)) #20443

# filter out NA values
x <- which(is.na(genes_2weeks$pvalue))
genes_2weeks <- genes_2weeks[-x,]
x <- which(is.na(genes_4weeks$pvalue))
genes_4weeks <- genes_4weeks[-x,]

# wilcox test
x_1 <- wilcox.test(x = genes_2weeks$pvalue[genes_2weeks$gene_names %in% secretome_sig$Gene.names],
                   y = genes_2weeks$pvalue[!genes_2weeks$gene_names %in% secretome_sig$Gene.names], alternative = "less", conf.int = TRUE)

y_1 <- wilcox.test(x = genes_4weeks$pvalue[genes_4weeks$gene_names %in% secretome_sig$Gene.names],
                   y = genes_4weeks$pvalue[!genes_4weeks$gene_names %in% secretome_sig$Gene.names], alternative = "less", conf.int = TRUE)

# 2 weeks timepoint filtered gene list
test <- which(!genes_2weeks$gene_names %in% secretome_sig$Gene.names)
test1 <- which(genes_2weeks$gene_names %in% secretome_sig$Gene.names)
genes_2_rna <- genes_2weeks[test,]
genes_2_sec <- genes_2weeks[test1,]

median(genes_2_rna$pvalue, na.rm = TRUE) # 0.05432413
median(genes_2_sec$pvalue, na.rm = TRUE) # 0.002349141

x <- as.data.frame(log2(genes_2_rna$baseMean + 1))
rownames(x) <- genes_2_rna$genes
y <- as.data.frame(log2(genes_2_sec$baseMean + 1))
rownames(y) <- genes_2_sec$genes
boxplot.default(x,y)
d1 <- density(genes_2_rna$pvalue)
d2 <- density(genes_2_sec$pvalue)
plot(d1)
plot(d2)

genes_2_rna$group <- as.factor("RNA-genes")
genes_2_sec$group <- as.factor("secretome-genes")
genes_2weeks <- rbind(genes_2_rna, genes_2_sec)
#genes_2weeks$logFC <- log2(genes_2weeks$baseMean)
plot_2weeks <- ggplot(genes_2weeks, aes(x=pvalue, color = group)) + geom_density() + xlim(0,1)

# 4 weeks timepoint filtered gene list
# -----------
test2 <- which(!genes_4weeks$gene_names %in% secretome_sig$Gene.names)
test3 <- which(genes_4weeks$gene_names %in% secretome_sig$Gene.names)
genes_4_rna <- genes_4weeks[test2,]
genes_4_sec <- genes_4weeks[test3,]

median(genes_4_rna$pvalue, na.rm = TRUE) #0.009765904
median(genes_4_sec$pvalue, na.rm = TRUE) # 2.680552e-05

z <- genes_4_rna$pvalue
q <-genes_4_sec$pvalue
boxplot.default(z,q)

d3 <- density(z)
d4 <- density(q)
plot(d3)
plot(d4)

genes_4_rna$group <- as.factor("RNA-genes")
genes_4_sec$group <- as.factor("secretome-genes")
genes_4weeks <- rbind(genes_4_rna, genes_4_sec)
genes_4weeks$logFC <- log2(genes_4weeks$baseMean)
plot_4weeks <- ggplot(genes_4weeks, aes(x=pvalue, color = group)) + geom_density() + xlim(0,1) + ylim(0, 30)
