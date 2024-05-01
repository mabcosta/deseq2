suppressMessages(library(EnhancedVolcano))
suppressMessages(library(tximeta))
suppressMessages(library(biomaRt))
suppressMessages(library(IHW))
suppressMessages(library(clusterProfiler))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(ReactomePA))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))
suppressMessages(library(DESeq2))
suppressMessages(library(apeglm))
suppressMessages(library(ggfortify))
suppressMessages(library(optparse))
suppressMessages(library(readr))
suppressMessages(library(dplyr))


## command line options
option_list <- list( 
        # specify output directory
        make_option(c('-o', '--outdir'), action='store', default=NA,
        type='character', help='output directory'),
        # choose first condition, also used as reference level
        make_option(c('-a', '--conditionA'), action='store', default=NA, type='character', help='first condition for comparison, also used as reference level for DESeq2'),
        # choose second condition for comparison
        make_option(c('-b', '--conditionB'), action='store', default=NA, type='character', help='second condition to be compared against condition_a'),
        # padj cut-off value - default is 0.05
        make_option(c('-c', '--padjCut'), action='store', default=0.05,type='double', help='padj cut-off value'),
        # volcano plot fold change (FC) cutoff - default is 0.5
        make_option(c('-x', '--VfcCut'), action='store',default=0.5, type='double', help='Fold change cutoff for Volcano'),
        # volcano plot pCutOff - default is 0.05
        make_option(c('-y', '--VpCut'), action='store', default=0.05, type='double', help='p-val cutoff for Volcano plot')
)

opt <- parse_args(OptionParser(option_list = option_list))

if (!is.na(opt$outdir)) {
  output_dir <- opt$outdir
  label <- basename(opt$outdir)
} else {
  print("error: no output directory specified")
}

get
# read in counts
countData <- read.csv("/home/mcosta/projects/deseq2/interger_counts_data_MIR184.csv", row.names = 1)

# read-in sample metadata (all)
colData <- read.csv("/home/mcosta/projects/deseq2/sample_info.csv", row.names = 1 )
# define conditions
coldata <- countData[condition %in% c(opt$conditionA,
                                   opt$conditionB)]

coldata <- countData[condition %in% c("scramble",
                                      "mir184")]



colData$condition <- factor(colData$condition)
coldata$condition <- factor(coldata$condition)

## DESeq2 standard workflow
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design= ~ condition)
keep <- rowSums(counts(dds) >= 5) >= 4
dds <- dds[keep,]
dds$condition <- relevel(dds$condition, ref = opt$conditionA)
dds <- DESeq(dds)
res <- results(dds, filterFun=ihw, alpha = opt$padjCut) # MLE

# shrinkage - generate MAP values for plotting and ranking only
coef <- resultsNames(dds)[2]
print(coef)
resLFC <- lfcShrink(dds, coef=coef, type="apeglm") # MAP

pdf(paste(output_dir, label, '_MA.pdf', sep =''))
plotMA(resLFC, ylim=c(-2,2))
dev.off()

# convert resLFC and res objects to data tables and Gencode IDs to Ensembl
resLFC$geneID <- gsub("\\..*","", rownames(resLFC))
res$geneID <- gsub("\\..*","", rownames(resLFC))
resLFC_dt <- as.data.table(resLFC)
res_dt <- as.data.table(res)

# merge MAP logfoldchange values to MLE table, rename the logfoldchange columns to .MLE and .MAP
map <- resLFC_dt[, .(log2FoldChange, geneID)] 
res1 <- merge(res_dt, map, by = 'geneID') # final DEGs table
setnames(res1, old=c("log2FoldChange.x","log2FoldChange.y"), new=c("log2FoldChange.MLE", "log2FoldChange.MAP"))

# add gene names to final DEGs table
genes <- res1$geneID
mart <- useEnsembl("genes", "hsapiens_gene_ensembl")
name <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol","description", "gene_biotype", "goslim_goa_description"),values=genes, mart= mart)
name.merged <- name %>%
  dplyr::group_by(ensembl_gene_id, hgnc_symbol, description, gene_biotype) %>%
  dplyr::summarise(goslim_goa_descriptions = paste(goslim_goa_description, collapse = ","))
res_named <- merge(res1,name.merged, by.x='geneID', by.y="ensembl_gene_id", sort = FALSE)

# create table of MAP values for plotting
genes.p <- resLFC_dt$geneID
name.p <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=genes.p, mart= mart)
plot_lfc <- merge(resLFC_dt,name.p, by.x='geneID', by.y="ensembl_gene_id", sort = FALSE) # plotting table

# sort and subset significant events, write to csv
resOrdered <- res_named[order(abs(res_named$log2FoldChange.MAP), decreasing=TRUE),]
fwrite(resOrdered, paste(output_dir, label, '_deseq2_results_ihw_no_filter.csv', sep =''))

resSig <- subset(resOrdered, padj < as.numeric(opt$padjCut))
fwrite(resSig, paste(output_dir, label, '_deseq2_results_ihw_padj_filtered.csv', sep =''))

print('results written')

# Data transformation and visualisation -  running times are shorter when using blind=FALSE and if the function DESeq() has already been run
vsd <- rlog(dds, blind = FALSE)

## loading genes
ntop <- 500
rv <- rowVars(assay(vsd))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
vsd_matrix <- t(assay(vsd)[select, ])

vsd_pca <- prcomp(vsd_matrix) # transpose so that PC is calculated for genes instead of sample
load <- as.data.frame(vsd_pca$rotation)
pc1 <- load[order(load$PC1, decreasing = TRUE),]
pc1 <- subset(pc1, select = PC1)
fwrite(pc1, paste(output_dir, label, '_pc1_loading_genes.csv', sep =''), row.names = TRUE)

pc2 <- load[order(load$PC2, decreasing = TRUE),]
pc2 <- subset(pc2, select = PC2)
fwrite(pc2, paste(output_dir, label, '_pc2_loading_genes.csv', sep =''), row.names = TRUE)

# PCA
title <- strsplit(coef, '_')[[1]]

svg(paste(output_dir, label, '_pca.svg', sep =''))
pcaData <- plotPCA(vsd, returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) + geom_point(size=1) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text_repel(label = coldata$names) + labs(title = paste(toupper(title[2]), 'vs', toupper(title[4]), 'with DESeq2')) + coord_fixed() + theme_bw()
dev.off()
print('PCA plot saved')

# # volcano plot
svg(paste(output_dir, label, '_volcano.svg', sep =''))
EnhancedVolcano(plot_lfc, lab = plot_lfc$hgnc_symbol, x = 'log2FoldChange', y = 'padj', ylab = expression(-log[10]~padj), legendPosition = 'none', drawConnectors = TRUE, max.overlaps = 7, title = paste(toupper(title[2]), 'vs', toupper(title[4])), pCutoff = as.numeric(opt$VpCut), FCcutoff = as.numeric(opt$VfcCut))
dev.off()
print('Volcano plot saved')

### Cluster Profiler ###
print('DESeq2 analysis complete')
print('Starting Cluster Profiler analysis')

# GO enrichment & GSEA
gene <- mapIds(org.Hs.eg.db, keys=resSig$hgnc_symbol, 
column="ENTREZID", keytype="SYMBOL", multiVals="first")

genelist <- resSig[,log2FoldChange.MAP]
names(genelist) = as.character(mapIds(org.Hs.eg.db, keys=resSig[,hgnc_symbol], column="ENTREZID", keytype="SYMBOL", multiVals="first"))
genelist <- genelist[!(names(genelist)) == 'NULL']
genelist <- sort(genelist, decreasing = TRUE)

ont <- list('MF', 'BP', 'CC')

for (i in ont) {
    print(i)
    go <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, ont = i, pAdjustMethod = 'BH', pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE) 
    gse <- gseGO(geneList = genelist, OrgDb = org.Hs.eg.db, ont = i, minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1, verbose = FALSE, nPermSimple = 10000)

    if (nrow(go) > 0) {
        fwrite(as.data.frame(go), paste(output_dir, label, '_go_enrichment_', i,'.csv', sep =''))
    } else {
        print('no GO enrichment')
    }

    if (nrow(gse) > 0) {
        fwrite(as.data.frame(gse), paste(output_dir, label, '_go_gsea_enrichment_', i,'.csv', sep =''))
    } else {
        print('no GSEA enrichment')
    }
}

rpa_path <- enrichPathway(gene = gene, organism = "human", pAdjustMethod = 'BH', pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)

if (nrow(rpa_path) > 0) {
        fwrite(as.data.frame(rpa_path), paste(output_dir, label, '_reactomePA_enrichment.csv', sep =''))
    } else {
        print('no reactome enrichment')
}


wiki_path <- enrichWP(gene = gene, pAdjustMethod = 'BH', pvalueCutoff = 1, qvalueCutoff = 1, organism = "Homo sapiens") 

if (nrow(wiki_path) > 0) {
        fwrite(as.data.frame(wiki_path), paste(output_dir, label, '_wikipaths_enrichment.csv', sep =''))
    } else {
        print('no wiki enrichment')
}


kegg <- enrichKEGG(gene = gene, organism = 'hsa', pvalueCutoff = 1, pAdjustMethod = 'BH', qvalueCutoff = 1)

mkegg <- enrichMKEGG(gene = gene, organism = 'hsa', pvalueCutoff = 1, pAdjustMethod = 'BH', qvalueCutoff = 1)

kegg_gsea <- gseKEGG(geneList = genelist, organism = 'hsa', minGSSize = 10, pvalueCutoff = 1, pAdjustMethod = "BH", verbose = FALSE, nPermSimple = 10000)

if (nrow(kegg) > 0) {
        fwrite(as.data.frame(setReadable(kegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")), paste(output_dir, label, '_kegg_enrichment.csv', sep =''))
    } else {
        print('no KEGG enrichment')
    }

if (nrow(as.data.frame(mkegg)) > 0) {
        fwrite(as.data.frame(setReadable(mkegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")), paste(output_dir, label, '_kegg_module_enrichment.csv', sep =''))
    } else {
        print('no KEGG module enrichment')
    }


if (nrow(kegg_gsea) > 0) {
    fwrite(as.data.frame(setReadable(kegg_gsea, OrgDb = org.Hs.eg.db, keyType="ENTREZID")), paste(output_dir, label, '_kegg_gsea_enrichment.csv', sep =''))
} else {
    print('no KEGG GSEA enrichment')
}