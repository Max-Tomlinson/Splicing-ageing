library(DESeq2)
library(ggplot2)
library(ggsci)
library(pheatmap)
library(RColorBrewer)
library(rrcov)
library(umap)

theme_new <-  theme_classic(base_size = 24) +
  theme(legend.position = "none", 
        text = element_text(family = "Helvetica"))
theme_set(new = theme_new)

countData <- read.csv(file =  "~/Desktop/Nanopore/Bambu/bambu2/counts_gene.txt", sep = "\t", row.names = 1)

colnames(x = countData) <- gsub(pattern = "^(.*?)_(.*?)_.*", replacement = "\\1_\\2", x = colnames(x = countData))

colData <- data.frame(colnames = colnames(x = countData))
colData$Time <- factor(x = ifelse(test = grepl(pattern = "^sample_0[0-9]|^sample_1[0-9]|^sample_2[0-5]", x = colData$colnames), yes = "Time_1", no = "Time_2"))

dds <- DESeqDataSetFromMatrix(countData = round(x = countData), colData = colData, design = as.formula(object = "~ Time"))
dds <- dds[rowSums(x = counts(object = dds) >= 10) >= 18,]
dds <- estimateSizeFactors(object = dds)

rld <- rlog(object = dds, blind = FALSE)
rld$Time <- gsub(pattern = "_", replacement = " ", x = rld$Time)
rld$pair <- rep(x = 1:18, 2)

pcaData <- plotPCA(object = rld, intgroup = c("Time", "pair"), returnData = TRUE)

umapData <- umap(d = t(x = assay(x = rld)))
umapData <- data.frame(UMAP1 = umapData$layout[, 1], UMAP2 = umapData$layout[, 2], Time = colData(rld)$Time, pair = colData(rld)$pair)
umapData$name <- rownames(x = umapData)

robustPCAData <- PcaGrid(x = as.matrix(x = t(x = assay(x = rld))))
robustPCAData <- data.frame(`Score distance` = robustPCAData$sd, `Score distance cut-off` = robustPCAData$cutoff.sd, Index = colnames(x = rld), Time = colData(rld)$Time, pair = colData(rld)$pair)

pca_plot <- ggplot(pcaData, mapping = aes(x = PC1, y = PC2, colour = Time, shape = factor(x = pair))) +
  geom_point(size = 5) +
  guides(colour = guide_legend(title = NULL), shape = "none") +
  scale_colour_futurama() +
  scale_shape_manual(values = 1:18)

umap_plot <- ggplot(umapData, mapping = aes(x = UMAP1, y = UMAP2, colour = Time, shape = factor(x = pair))) +
  geom_point(size = 5) +
  guides(colour = guide_legend(title = NULL), shape = "none") +
  scale_colour_futurama() +
  scale_shape_manual(values = 1:18)

rpc_plot <- ggplot(robustPCAData, mapping = aes(x = factor(x = Index), y = Score.distance, colour = Time, shape = factor(x = pair))) +
  geom_hline(yintercept = mean(x = robustPCAData$Score.distance.cut.off), colour = "red") +
  geom_point(size = 5) +
  guides(colour = guide_legend(title = NULL), shape = "none") +
  labs(title = "Robust principal component analysis", x = "Samples", y = "Score distance") +
  scale_colour_futurama() +
  scale_shape_manual(values = 1:18) +
  theme(axis.text.x = element_blank())

genes <- order(rowMeans(x = counts(object = dds, normalized = TRUE)), decreasing = TRUE)[1:500]

annotation_col <- as.data.frame(x = colData(x = rld)["Time"])
annotation_col$Time <- sub(pattern = "Time ", replacement = "", x = annotation_col$Time)
annotation_colors <- list(Time = c("1" = "orange", "2" = "red"))
  
pheatmap(mat = assay(x = rld)[genes, ], color = colorRampPalette(colors = brewer.pal(n = 9, name = "YlOrRd"))(100), annotation_colors = annotation_colors, annotation_legend = FALSE, show_colnames = FALSE, show_rownames = FALSE, annotation_col = annotation_col, main = "Heatmap of the top 500 most variable genes", fontsize = 22, filename = "~/Desktop/Nanopore/Analysis/PCA/Pheatmap.pdf", width = 11, height = 8.5)

ggsave(filename = "~/Desktop/Nanopore/Analysis/PCA/PCA.pdf", plot = pca_plot, width = 11, height = 8.5)
ggsave(filename = "~/Desktop/Nanopore/Analysis/PCA/UMAP.pdf", plot = umap_plot, width = 11, height = 8.5)
ggsave(filename = "~/Desktop/Nanopore/Analysis/PCA/RPC.pdf", plot = rpc_plot, width = 11, height = 8.5)
