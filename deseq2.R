setwd("/home/mcosta/MarcosCosta/projects/MIR184/DESeq2")

library(DESeq2)
library(tidyverse)
library(apeglm)

# Prepare counts data

# read in counts
countData <- read.csv("interger_counts_data_MIR184.csv", row.names = 1)
head(countData)
ncol(countData)
# read in sample info

colData <- read.csv("sample_info.csv", row.names = 1 )
head(colData,10)
nrow(colData)
# Making sure row names in colData match and are in the same order as counts_data
all(colnames(countData) %in% rownames(colData))
all(colnames(countData) == rownames(colData))

# Example code to convert a character variable to a factor
colData$condition <- factor(colData$condition)

# Set up DESeq2 matrix
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design= ~ condition )

# Pre-filtering to remove low gene counts and NAs - improves graph visualisation
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,] # Dimension reduced from 63241x24 to 18923x24

# Run DESeq2 - takes less than a minute
dds <- DESeq(dds)

# Analysis contrast between conditions
resultsNames(dds) # lists the coefficients
res <- results(dds, name="condition_mir184_vs_anti_mir184")
res

# Shrink log fold changes using apeglm 
resLFC <- lfcShrink(dds, coef="condition_mir184_vs_anti_mir184", type="apeglm")
resLFC

# Order results by p-value
resOrdered <- res[order(res$pvalue),]
summary(res)
sum(res$padj < 0.1, na.rm=TRUE)

res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res$padj < 0.05, na.rm=TRUE)


# Independent hypothesis weighting
library("IHW")
resIHW <- results(dds, filterFun=ihw)

summary(resIHW)
sum(resIHW$padj < 0.1, na.rm=TRUE)
metadata(resIHW)$ihwResult

# Visualisation 
plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))

plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

# Getting heatmaps and PCAs
# Variance
ntd <- normTransform(dds)
vsd <- vst(dds, blind=FALSE)
head(assay(vsd), 3)
library("vsn")
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))

sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

plotPCA(vsd, intgroup=c("condition"))


# HTML report
library("ReportingTools")

des2Report <- HTMLReport(shortName = "RNAseq_analysis_with_DESeq2",
                         title = "RNA-seq using DESeq2",
                         reportDirectory = "./reports")

publish(dds,des2Report, pvalueCutoff=0.05, annotation.db="org.Mm.eg.db",
        factor = colData(dds)$conditions, reportDir="./reports")

finish(des2Report)
