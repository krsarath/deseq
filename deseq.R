library(BiocManager)
library(DESeq2)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(plotly)
library(RColorBrewer)
library(EnhancedVolcano)
library(ggVennDiagram)
library(tidyverse)

setwd('/home/sarath/R/rwd/deseqnew')

count_data <- read.csv('matrixcounts.csv',header = TRUE,row.names = 1)
View(count_data)
colnames(count_data)
head(count_data)

sample_info <- read.csv('sample.csv',header = TRUE)
colnames(sample_info)
head(sample_info)

sample_info$Treatment <- factor(sample_info$Treatment)

countdata <- count_data[which(rowMeans(!is.na(count_data)) > 0.5), ]
head(countdata)
countdata <- countdata %>% mutate(across(where(is.numeric), ~ replace(., is.na(.), median(., na.rm = TRUE))))
head(countdata)
write.csv(countdata,'COUNTDATA.csv')

dds <- DESeqDataSetFromMatrix(countData = round(countdata), colData = sample_info, design = ~ 0 + Treatment)
dds$Treatment <- factor(dds$Treatment, levels = c('control','treated'))
keep <- rowSums(counts(dds) > 10) >= min(table(sample_info$Treatment))
dds <- dds[keep,]
dds <- DESeq(dds, fitType = 'local')
deseq_result <- results(dds)
head(deseq_result)
View(deseq_result)

deseq_result <- as.data.frame(deseq_result)
class(deseq_result)
head(deseq_result)

deseq_result_ordered <- deseq_result[order(deseq_result$pvalue),]
head(deseq_result_ordered)
deseq_result['TRINITY_DN27872_c0_g1_i2',]
write.csv(deseq_result,'deresult.csv')

filtered <- deseq_result %>% filter(deseq_result$padj < 0.05)
filtered <- filtered %>% filter(abs(filtered$log2FoldChange) > 1)
dim(deseq_result)
dim(filtered)
write.csv(filtered,'deresfiltered.csv')

plotDispEsts(dds)
vsd <- vst(dds,blind = F, fitType = 'local')
plotPCA(vsd,intgroup= 'Treatment')
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
colnames(sampleDistMatrix)
write.csv(sampleDistMatrix, 'distancematrix.csv')

colors <- colorRampPalette(rev(brewer.pal(9,"Blues")))(255)
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists, col=colors)

top_hits <- deseq_result[order(deseq_result$padj),][1:10,]
top_hits <- row.names(top_hits)
top_hits

rld <-rlog(dds,blind = F)
pheatmap(assay(rld)[top_hits,],cluster_rows = F,show_rownames = T,cluster_cols = F)
pheatmap(assay(rld)[top_hits,])
annot_info <- as.data.frame(colData(dds))
heatmap_data <- assay(rld)[top_hits, ]
summary(heatmap_data)
write.csv(summary(heatmap_data), 'summaryheatmapdata.csv')

plotMA(dds,ylim=c(-12,12))

resultsNames(dds)

resLFC <- lfcShrink(dds,contrast = c("Treatment","control","treated"), type='ashr')
write.csv(resLFC,"reslfc.csv")
plotMA(resLFC,ylim=c(-10,10))
head(resLFC)
resordered <- resLFC[order(resLFC$pvalue),]
head(resordered)
summary(resordered)
write.csv(resordered,"resorderd.csv")
sum(resLFC$padj < 0.1, na.rm=TRUE)
NAresfilt <- results(dds, alpha=0.05)
summary(NAresfilt)

hist(NAresfilt$padj)
resLFC <- as.data.frame(resLFC)
resLFC$diffexpressed <- NA
resLFC$diffexpressed[resLFC$log2FoldChange > 1 & resLFC$padj < 0.05] <- "UP"
resLFC$diffexpressed[resLFC$log2FoldChange < -1 & resLFC$padj < 0.05] <- "DOWN"
resLFC$delabel <- NA

colnames(resLFC)
head(resLFC)

ggplot(data = resLFC,aes(x = log2FoldChange,y = -log10(pvalue),col=diffexpressed,label=delabel))+
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2)+
  geom_text_repel()+
  scale_color_manual(values=c('blue','red'),labels = c("Downregulated", "Not significant", "Upregulated"))+
  theme(text = element_text(size = 20))

g1 <- ggVennDiagram(deseq_result, category.names = c('baseMean', 'log2FoldChange', 'lfcSE', 'pvalue', 'padj'), show_intersect = TRUE, label = "both", label_alpha = 0)
g1 + scale_x_continuous(expand = expansion(mult = .2))
g1 + scale_fill_distiller(palette = "Greens")
g1

