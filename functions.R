# Functions that are used for the Ska2 project from Jakob Hartmann

# Converting ensembl names to external gene names
biomart_ensembl_to_gene <- function(gene){
  ensembl=useMart("ensembl")
  ensembl = useDataset("mmusculus_gene_ensembl",mart=ensembl)
  names <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                 filters = 'ensembl_gene_id', 
                 values = gene, 
                 mart = ensembl)
  return(names)
}

# This function not just replaces SVD with PCA, but also performs an anova(lm())
# test for both numeric and factor variables, instead of an lm and kruskal test,
# as done by the original function. It also adds percentage of variance 
# explained for each PC to the axis label, and the F value (a measure of effect 
# size) to each box. Besides the centering of the beta values as done by the 
# original function, this function also scales and centers the numeric variables
# of pd.

library(genefilter)
library(ggplot2)
library(ggrepel)
library("factoextra")

custom.PCA <- function (beta, rgSet = NULL, pd, RGEffect = FALSE, 
                        PDFplot = TRUE, Rplot = TRUE, 
                        resultsDir = "./PCAimages/",
                        cex.axis = 1,
                        adjust.p = FALSE,
                        return.p = FALSE,
                        npca = NULL,
                        plot.title) 
{
  message("[===========================]")
  message("[<<<<< CUSTOM.PCA START >>>>>]")
  message("-----------------------------")
  GenPlot <- function(thdens.o, estdens.o, evalues.v) {
    minx <- min(min(thdens.o$lambda), min(evalues.v))
    maxx <- max(max(thdens.o$lambda), max(evalues.v))
    miny <- min(min(thdens.o$dens), min(estdens.o$y))
    maxy <- max(max(thdens.o$dens), max(estdens.o$y))
  }
  EstDimRMTv2 <- function(data.m) {
    M <- data.m
    for (c in 1:ncol(M)) M[, c] <- (data.m[, c] - mean(data.m[, c], na.rm = TRUE))/sqrt(var(data.m[, c], na.rm = TRUE))
    sigma2 <- var(as.vector(M), na.rm = TRUE)
    Q <- nrow(data.m)/ncol(data.m)
    thdens.o <- thdens(Q, sigma2, ncol(data.m))
    C <- 1/nrow(M) * t(M) %*% M
    eigen.o <- eigen(C, symmetric = TRUE)
    estdens.o <- density(eigen.o$values, from = min(eigen.o$values), 
                         to = max(eigen.o$values), cut = 0)
    GenPlot(thdens.o, estdens.o, eigen.o$values)
    intdim <- length(which(eigen.o$values > thdens.o$max))
    return(list(cor = C, dim = intdim, estdens = estdens.o, 
                thdens = thdens.o))
  }
  thdens <- function(Q, sigma2, ns) {
    lambdaMAX <- sigma2 * (1 + 1/Q + 2 * sqrt(1/Q))
    lambdaMIN <- sigma2 * (1 + 1/Q - 2 * sqrt(1/Q))
    delta <- lambdaMAX - lambdaMIN
    roundN <- 3
    step <- round(delta/ns, roundN)
    while (step == 0) {
      roundN <- roundN + 1
      step <- round(delta/ns, roundN)
    }
    lambda.v <- seq(lambdaMIN, lambdaMAX, by = step)
    dens.v <- vector()
    ii <- 1
    for (i in lambda.v) {
      dens.v[ii] <- (Q/(2 * pi * sigma2)) * sqrt((lambdaMAX - i) * (i - lambdaMIN))/i
      ii <- ii + 1
    }
    return(list(min = lambdaMIN, max = lambdaMAX, step = step, 
                lambda = lambda.v, dens = dens.v))
  }
  runpca <- function(featuresbysamples, npca) { # modified function from DEGreport
    pca.res <- prcomp(t(featuresbysamples), center = FALSE, retx = TRUE)
    pc.var <- pca.res$sdev^2L
    pve <- 100L * (pc.var/sum(pc.var))
    #npca <- max(1L, length(which(pve >= min_pve_pct_pc)))
    samplepcvals <- pca.res$x[, 1L:npca, drop = FALSE]
    list(samplepcvals = samplepcvals, pve = pve)
  }
  if ((!file.exists(resultsDir)) & PDFplot) 
    dir.create(resultsDir)
  if (PDFplot) 
    message("custom.PCA Results will be saved in ", resultsDir, " .\n")
  if (length(which(is.na(beta))) > 0) 
    message(length(which(is.na(beta))), 
            " NAs are detected in your beta value data set, which may cause errors in the PCA. You may want to impute NAs first.")
  message("[PCA analysis will proceed with ", dim(beta)[1], " probes and ", 
          dim(beta)[2], " samples.]\n")
  message("\n[ custom.PCA() will only check the dimensions between data and pd, instead of checking if sample names are correctly matched, thus please make sure your pd file is in accord with your data sets (beta values and rgSet).]\n")
  if (is.null(pd) | class(pd) == "list") 
    stop("pd parameter in data frame or matrix format is necessary and must contain at least two factors. If your pd is a list, please change its format.")
  if (class(pd) == "matrix") 
    pd <- as.data.frame(pd)
  PhenoTypes.lv_tmp <- pd[, !colnames(pd) %in% c("filenames", "Basename") & apply(pd, 2, function(x) length(unique(x))) != 1]
  PhenoTypes.lv <- PhenoTypes.lv_tmp
  if (!is.null(rownames(pd))) 
    rownames(PhenoTypes.lv) <- rownames(pd)
  if (ncol(PhenoTypes.lv) >= 2) {
    message("<< The following factors in your pd will be analyzed: >>")
    sapply(colnames(PhenoTypes.lv_tmp), 
           function(x) message("<", x, ">(", class(PhenoTypes.lv[[x]]), "):", paste(unique(PhenoTypes.lv_tmp[, x]), collapse = ", ")))
    message("[custom.PCA has automatically selected ALL factors that contain at least two different values from pd. If you don't want to analyze some of them, please remove them from pd.]")
  }
  else {
    stop("There are no factors with at least two values to be analyzed. Maybe your factors contain only one value?")
  }
  if (ncol(pd) > ncol(PhenoTypes.lv)) {
    message("\n<< The following factors in pd will not be analyzed: >>")
    sapply(setdiff(colnames(pd), colnames(PhenoTypes.lv)), 
           function(x) message("<", x, ">"))
    message("[Factors are ignored because they only indicate name or project, or they contain only one value across all samples.]")
  }
  if (RGEffect == TRUE & is.null(rgSet)) 
    message("If you want to check the effects of control probes, you must provide an rgSet. Now custom.SVD can only analyze factors in pd.")
  if (!is.null(rgSet) & RGEffect) {
    if (rgSet@annotation[1] == "IlluminaHumanMethylation450k") 
      data(ControlProbes450K)
    else data(ControlProbesEPIC)
    dataC2.m <- as.data.frame(log2(apply(ControlProbes, 1, 
                                         function(x) if (x[3] == "Grn") getGreen(rgSet)[x[2],] else getRed(rgSet)[x[2], ])))
    PhenoTypes.lv <- cbind(PhenoTypes.lv, dataC2.m)
    message("\n<< The following rgSet information has been included: >>")
    sapply(colnames(dataC2.m), function(x) message("<", x, ">"))
  }
  if (nrow(PhenoTypes.lv) == ncol(beta)) 
    message("\n<< PhenoTypes.lv generated successfully. >>")
  else stop("Dimension of your pd file (and rgSet information) is not equal to your beta value matrix.")
  beta <- as.matrix(beta)
  tmp.m <- beta - rowMeans(beta, na.rm = TRUE)
  if(is.null(npca)){
    npca <- EstDimRMTv2(tmp.m)$dim
  }
  pca.o <- runpca(tmp.m, npca)
  if (ncol(pca.o$samplepcvals) > 20) 
    topPCA <- 20
  else topPCA <- ncol(pca.o$samplepcvals)
  pcaPV.m <- matrix(nrow = topPCA, ncol = ncol(PhenoTypes.lv))
  pcaFV.m <- matrix(nrow = topPCA, ncol = ncol(PhenoTypes.lv))
  colnames(pcaPV.m) <- colnames(pcaFV.m) <- colnames(PhenoTypes.lv)
  # for (c in 1:topPCA) for (f in 1:ncol(PhenoTypes.lv)) if (class(PhenoTypes.lv[, f]) != "numeric") 
  #     pcaPV.m[c, f] <- kruskal.test(pca.o$samplepcvals[, c] ~ as.factor(PhenoTypes.lv[[f]]))$p.value
  # else pcaPV.m[c, f] <- summary(lm(pca.o$samplepcvals[, c] ~ PhenoTypes.lv[[f]]))$coeff[2, 4]
  for (c in 1:topPCA) for (f in 1:ncol(PhenoTypes.lv)) if (!any(class(PhenoTypes.lv[, f]) %in% c("numeric", "integer"))) {
    pcaPV.m[c, f] <- anova(lm(pca.o$samplepcvals[, c] ~ as.factor(PhenoTypes.lv[[f]])))$`Pr(>F)`[1]
    pcaFV.m[c, f] <- round(anova(lm(pca.o$samplepcvals[, c] ~ as.factor(PhenoTypes.lv[[f]])))$`F value`[1], 0)
  } else {
    pcaPV.m[c, f] <- anova(lm(pca.o$samplepcvals[, c] ~ scale(PhenoTypes.lv[[f]])))$`Pr(>F)`[1]
    pcaFV.m[c, f] <- round(anova(lm(pca.o$samplepcvals[, c] ~ scale(PhenoTypes.lv[[f]])))$`F value`[1], 0)
  }
  if(adjust.p){
    pcaPV.m <- matrix(p.adjust(pcaPV.m, method = "bonferroni"), nrow = nrow(pcaPV.m), 
                      dimnames = list(rownames(pcaPV.m), colnames(pcaPV.m)))
  }
  pcaFV.m[pcaFV.m > 999] <- ">999" # prevents printing huge numbers that do not fit in the cells
  message("<< Calculated PC matrix successfully. >>")
  myPalette <- c("darkred", "red", "orange", "pink", "white")
  breaks.v <- c(-200, -10, -5, -2, log10(0.05), 0)
  if (Rplot) {
    par(mar = c(5, 15, 2, 1))
    image(x = 1:nrow(pcaPV.m), y = 1:ncol(pcaPV.m), z = log10(pcaPV.m + .Machine$double.eps), # added smallest fraction to prevent generation of infinite values 
          col = myPalette, breaks = breaks.v, xlab = "", ylab = "", 
          axes = FALSE, main = plot.title)
    text(x = row(pcaFV.m), y = col(pcaFV.m), pcaFV.m, cex = 0.5)
    axis(1, at = 1:nrow(pcaPV.m), labels = paste0("PC-", 1:nrow(pcaPV.m), 
                                                  " (", round(pca.o$pve[1:topPCA], 2), 
                                                  "%)"), 
         las = 2, cex.axis = 0.6)
    suppressWarnings(axis(2, at = 1:ncol(pcaPV.m), 
                          labels = colnames(pcaPV.m), las = 2, 
                          cex.axis = cex.axis))
    legend(x = -(topPCA/2.2), y = 3, 
           legend = c(expression("p < 1x" ~ 10^{-10}), expression("p < 1x" ~ 10^{-5}), "p < 0.01", "p < 0.05", "p > 0.05"), 
           fill = c("darkred", "red", "orange", "pink", "white"), 
           par("usr")[2], 
           par("usr")[4], xpd = NA,
           bty = "n",
           cex = 0.7)
  }
  if (PDFplot) {
    pdf(paste(resultsDir, "PCAsummary.pdf", sep = ""), width = 10, 
        height = 8)
    par(mar = c(5, 15, 2, 1), xpd = TRUE)
    image(x = 1:nrow(pcaPV.m), y = 1:ncol(pcaPV.m), z = log10(pcaPV.m), 
          col = myPalette, breaks = breaks.v, xlab = "", ylab = "", 
          axes = FALSE, main = plot.title)
    text(x = row(pcaFV.m), y = col(pcaFV.m), pcaFV.m, cex = 0.5)
    axis(1, at = 1:nrow(pcaPV.m), labels = paste0("PC-", 1:nrow(pcaPV.m), 
                                                  " (", round(pca.o$pve[1:topPCA], 2), 
                                                  "%)"), 
         las = 2, cex.axis = 0.6)
    suppressWarnings(axis(2, at = 1:ncol(pcaPV.m), 
                          labels = colnames(pcaPV.m), las = 2, 
                          cex.axis = cex.axis))
    legend(x = -(topPCA/2.2), y = 3, 
           legend = c(expression("p < 1x" ~ 10^{-10}), expression("p < 1x" ~ 10^{-5}), "p < 0.01", "p < 0.05", "p > 0.05"), 
           fill = c("darkred", "red", "orange", "pink", "white"), 
           par("usr")[2], 
           par("usr")[4], xpd = NA,
           bty = "n",
           cex = 0.7)
    dev.off()
  }
  message("<< Plotted PCA matrix successfully. >>")
  message("[<<<<<< custom.PCA END >>>>>>]")
  message("[===========================]")
  if(return.p){
    pcaPV.m <- matrix(apply(pcaPV.m, 2, function(x) sprintf("%.2e", x)), nrow = topPCA, dimnames = list(paste0("PC", 1:topPCA), colnames(pcaPV.m)))
    return(t(pcaPV.m)[ncol(pcaPV.m):1,, drop=FALSE])
  } 
}

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


plotPCA.deseq2_modified <- function (object, intgroup = "condition", 
                                     ntop = 500, returnData = FALSE) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  
  ## Select the PCAs and percentVar that you like instead of 1 and 2
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3],
                  PC4 = pca$x[, 4], PC5 = pca$x[, 5], PC6 = pca$x[, 6],
                  PC7 = pca$x[, 7], PC8 = pca$x[, 8], PC9 = pca$x[, 9],
                  PC10 = pca$x[, 10], PC11 = pca$x[, 11], 
                  PC12 = pca$x[, 12],
                  PC13 = pca$x[, 13], PC14 = pca$x[, 14], 
                  PC15 = pca$x[, 15],
                  group = group, 
                  intgroup.df, name=colnames(object))
  
  
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:15]
    return(list(d,pca))
  }
  ggplot(data=d, aes_string(x="PC1", y="PC2", color="group"))
  + geom_point(size=3) + 
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
    ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
    coord_fixed()
  
}
