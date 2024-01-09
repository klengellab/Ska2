# Wilcoxon analysis Ska2 significant genes to compare against the secretome

# reference: https://stat.ethz.ch/R-manual/R-devel/library/stats/html/wilcox.test.html

#Libraries
library(ggplot2)
library(biomaRt)


# Data load
secretome <- read.csv("/PHShome/je637/RNAseq/RNAseq_Ska2/data/secretome_JH_2020_dec.csv", header = TRUE, sep = ",")
genes_2weeks <- read.csv("/PHShome/je637/RNAseq/RNAseq_Ska2/output/significant_genes_2weeks.csv")
genes_4weeks <- read.csv("/PHShome/je637/RNAseq/RNAseq_Ska2/output/significant_genes_4weeks.csv")

# Filter out significant genes on the secretome
secretome_sig <- secretome[,c(7,2)]
x <- which(secretome_sig$P.value < 0.05)
secretome_sig <- secretome_sig[x,]
# In total there are 446 genes significant different (p.value < 0.05) in the secretome data

# biomart to convert ensmbl id to external gene names
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

gene_names_2weeks <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                           filters = "ensembl_gene_id",
                           values = genes_2weeks$genes,
                           mart = ensembl)

gene_names_4weeks <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                           filters = "ensembl_gene_id",
                           values = genes_4weeks$genes,
                           mart = ensembl)

x <- match(gene_names_2weeks$ensembl_gene_id, genes_2weeks$genes)
y <- match(gene_names_4weeks$ensembl_gene_id, genes_4weeks$genes)

genes_2weeks <- genes_2weeks[x,]
genes_4weeks <- genes_4weeks[y,]

genes_2weeks$gene_names <- gene_names_2weeks$external_gene_name
genes_4weeks$gene_names <- gene_names_4weeks$external_gene_name

# Wilcoxon rank sum test

length(which(secretome_sig$Gene.names %in% genes_2weeks$gene_names)) #121
length(which(secretome_sig$Gene.names %in% genes_4weeks$gene_names)) #184

length(genes_2weeks$padj) == length(unique(genes_2weeks$padj))
length(unique(genes_2weeks$padj)) #3339
# There are quit some genes with the same padj therefore I am using the nominal p.values
length(unique(genes_2weeks$pvalue)) #3339

length(genes_4weeks$padj) == length(unique(genes_4weeks$padj))
length(unique(genes_4weeks$padj)) #5602
length(unique(genes_4weeks$pvalue)) #5602


# wilcox test
x_1 <- wilcox.test(x = genes_2weeks$padj[genes_2weeks$gene_names %in% secretome_sig$Gene.names],
                   y = genes_2weeks$padj[!genes_2weeks$gene_names %in% secretome_sig$Gene.names], alternative = "less", conf.int = TRUE)

y_1 <- wilcox.test(x = genes_4weeks$padj[genes_4weeks$gene_names %in% secretome_sig$Gene.names],
                   y = genes_4weeks$padj[!genes_4weeks$gene_names %in% secretome_sig$Gene.names], alternative = "less", conf.int = TRUE)
# 2 weeks timepoint filtered gene list but the basemeans are not normalized

test <- which(!genes_2weeks$gene_names %in% secretome_sig$Gene.names)
test1 <- which(genes_2weeks$gene_names %in% secretome_sig$Gene.names)
genes_2_rna <- genes_2weeks[test,]
genes_2_sec <- genes_2weeks[test1,]

median(genes_2_rna$pvalue, na.rm = TRUE) # 2.501842e-10
median(genes_2_sec$pvalue, na.rm = TRUE) # 7.633584e-12

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
plot_2weeks <- ggplot(genes_2weeks, aes(x=padj, color = group)) + geom_density() + xlim(0,1)
pdf(file = "/PHShome/je637/RNAseq/RNAseq_Ska2/output/density_wilcox_2.pdf")
plot_2weeks
dev.off()

# 4 weeks timepoint filtered gene list but the basemeans are not normalized
# -----------
test2 <- which(!genes_4weeks$gene_names %in% secretome_sig$Gene.names)
test3 <- which(genes_4weeks$gene_names %in% secretome_sig$Gene.names)
genes_4_rna <- genes_4weeks[test2,]
genes_4_sec <- genes_4weeks[test3,]

median(genes_4_rna$pvalue, na.rm = TRUE) #1.475277e-11
median(genes_4_sec$pvalue, na.rm = TRUE) #9.606171e-14

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
plot_4weeks <- ggplot(genes_4weeks, aes(x=padj, color = group)) + geom_density() + xlim(0,1) +theme_classic()
pdf(file = "/PHShome/je637/RNAseq/RNAseq_Ska2/output/density_wilcox_4.pdf")
plot_4weeks
dev.off()
