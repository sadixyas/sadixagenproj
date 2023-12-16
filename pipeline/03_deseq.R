#!/usr/bin/env Rscript
library(DESeq2)
library(tximport)
library(dplyr)
library(ggplot2)
library(magrittr)
library(Biobase)
library(pheatmap)
library(RColorBrewer)

samples <- read.table("samples.csv",header=TRUE,sep=",")
samples
samples$Name = sprintf("%s.r%s",samples$CONDITION,samples$REPLICATE)
samples$Name
files <- file.path("/bigdata/stajichlab/sadikshs/gen220/genproject/results/kallisto",samples$Name,"abundance.tsv")
names(files) <- samples$Name
#if (all(file.exists(files))) {
  
  txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)
  head(txi.kallisto$counts)
  head(txi.kallisto$abundance)
  #colnames(txi.kallisto$counts) <- samples$Name
  #colnames(txi.kallisto$abundance) <- samples$Name
  write.csv(txi.kallisto$abundance,"kallisto.TPM.csv")
  write.csv(txi.kallisto$counts,"kallisto.counts.csv")
  
  # DEseq2 analyses
  #geno = samples$CONDITION
  #treatment = samples$CONDITION
  #sampleTable <- data.frame(condition=treatment)
  
  #sampleTable <- data.frame(condition=treatment)
  samples$CONDITION <- factor(samples$CONDITION)
 # sampleTable$genotype <- factor(sampleTable$genotype)
  rownames(samples) = samples$Name
  
  dds <- DESeqDataSetFromTximport(txi.kallisto,samples, ~ CONDITION)
  
  #dds <- DESeqDataSetFromTximport(txi.kallisto,sampleTable, ~ condition)
  nrow(dds)
  dds <- dds[ rowSums(counts(dds)) > 1, ] #filter for null data
  nrow(dds)
  
  dds <- estimateSizeFactors(dds)
  vsd <- vst(dds, blind = FALSE) #normalization for plotting 
  
  #rld <- rlog(dds, blind = TRUE)
  #vsd <- varianceStabilizingTransformation(dds, blind = FALSE, fitType = "parametric")
  
  head(assay(vsd), 3)
  
  df <- bind_rows(
    as_tibble(log2(counts(dds, normalized=FALSE)[, 1:2]+1)) %>%
      mutate(transformation = "log2(x + 1)"), #part of tidyverse for reating a new column
    #as_tibble(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))
    as_tibble(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"))
  
  colnames(df)[1:2] <- c("x", "y")
  
  #for the first figure
  ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
    coord_fixed() + facet_grid( . ~ transformation)
  
  select <- order(rowMeans(counts(dds,normalized=FALSE)),
                  decreasing=TRUE)[1:30]
  df <- as.data.frame(colData(dds)[,c("CONDITION")])
  #for the second figure
  pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,cluster_cols=FALSE,annotation_col=samples)
  
  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(vsd$CONDITION, sep="-")
  #rownames(sampleDistMatrix) <- vsd$condition
  colnames(sampleDistMatrix) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  #for the third figure
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
  
  pcaData <- plotPCA(vsd, intgroup=c("CONDITION"), returnData=TRUE)
  #pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  #for the fourth figure
  ggplot(pcaData, aes(PC1, PC2, color=CONDITION, shape=name)) +
    #ggplot(pcaData, aes(PC1, PC2, color=condition)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed()
