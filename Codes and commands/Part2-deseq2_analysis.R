library(DESeq2)
library(tidyverse)
library(RColorBrewer)


counts_data <- read.table("gene_counts.tsv", header = TRUE, sep = "\t", row.names = 1)
colData <- read.table("sample_metadata.tsv", header = TRUE, sep = "\t", row.names = 1)
colData <- colData[colnames(counts_data), , drop = FALSE]
head(colData)
all(colnames(counts_data) == rownames(colData))
colData$condition <- factor(colData$condition)
str(colData)
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = colData, design = ~ condition)
dds
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds
dds$condition <- relevel(dds$condition, ref = "normal")
levels(dds$condition)
dds <- DESeq(dds)
res <- results(dds)
res <- res[order(res$padj), ]
summary(res)
res0.01 <- results(dds, alpha = 0.01)
res0.01 <- res0.01[order(res0.01$padj), ]
summary(res0.01)
plotMA(res)


# Top 10 UP-regulated genes
up_genes <- res0.01[which(res0.01$padj < 0.01 & res0.01$log2FoldChange > 0), ]
up_genes <- up_genes[order(-up_genes$log2FoldChange), ]
head(up_genes, 10)

# Top 10 DOWN-regulated genes
down_genes <- res0.01[which(res0.01$padj < 0.01 & res0.01$log2FoldChange < 0), ]
down_genes <- down_genes[order(down_genes$log2FoldChange), ]
head(down_genes, 10)
# topmost up-regulated gene
res0.01[which.max(res0.01$log2FoldChange), ]
# topmost down-regulated gene
res0.01[which.min(res0.01$log2FoldChange), ]
plotCounts(dds, gene = rownames(res)[which.min(res$padj)], intgroup = "condition")

# Plot most up-regulated gene
plotCounts(dds, gene = rownames(res0.01)[which.max(res0.01$log2FoldChange)], intgroup = "condition")
# Plot most down-regulated gene
plotCounts(dds, gene = rownames(res0.01)[which.min(res0.01$log2FoldChange)], intgroup = "condition")

mcols(res)$description
table(res$padj < 0.01)
res <- res[order(res$padj), ]
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized = TRUE)), 
                 by = "row.names", sort = FALSE)
names(resdata)[1] <- "Gene"
head(resdata)

write.table(resdata, file = "diffexpr-results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
hist(res$padj, breaks = 50, col = "grey", main = "Histogram of adjusted p-values")

# histogram log2FoldChange
hist(res$log2FoldChange, breaks = 50, col = "grey", main = "Histogram of log2 Fold Change")
rld <- rlogTransformation(dds)
head(assay(rld))



mycols <- brewer.pal(8, "Dark2")[1:length(unique(colData$condition))]

# Sample distance heatmap
sampleDists <- as.matrix(dist(t(assay(rld))))
library(gplots)
heatmap.2(as.matrix(sampleDists), key = FALSE, trace = "none", 
          col = colorpanel(100, "black", "white"),
          ColSideColors = mycols[colData$condition], 
          RowSideColors = mycols[colData$condition],
          margin = c(10, 10), main = "Sample Distance Matrix")


plotPCA(rld, intgroup = c("condition"))


sig_genes <- res0.01[which(res0.01$padj < 0.01), ]
sig_gene_ids <- rownames(sig_genes)
length(sig_gene_ids)

write.table(data.frame(gene_id = sig_gene_ids, 
                       log2FC = sig_genes$log2FoldChange,
                       padj = sig_genes$padj),
            file = "significant_genes.tsv", 
            sep = "\t", quote = FALSE, row.names = FALSE)
up_gene_ids <- rownames(res0.01[which(res0.01$padj < 0.01 & res0.01$log2FoldChange > 0), ])
down_gene_ids <- rownames(res0.01[which(res0.01$padj < 0.01 & res0.01$log2FoldChange < 0), ])

writeLines(up_gene_ids, "upregulated_genes.txt")
writeLines(down_gene_ids, "downregulated_genes.txt")
