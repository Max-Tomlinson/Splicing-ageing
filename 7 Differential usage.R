library(DEXSeq)
library(DRIMSeq)
library(ggplot2)
library(gprofiler2)
library(stageR)

dxr <- readRDS(file = "~/Desktop/Nanopore/Analysis/Differential usage/dxr.rds")

#DEXSeqHTML(object = dxr, path = "~/Desktop/Nanopore/Analysis/Differential usage", FDR = 0.1, colour = c("orange", "red"))

pdf(file = "~/Desktop/Nanopore/Analysis/Differential usage/DEXSeq_MBNL3.pdf", width = 11.5, height = 8)

plotDEXSeq(object = dxr, geneID = dxr[order(dxr$pvalue), ]$groupID[1], expression = FALSE, splicing = TRUE, displayTranscripts = FALSE, legend = TRUE, color = c("#4040FF", "#FF4162"), cex.axis = 1.2, cex = 2, lwd = 3)

dev.off()

pdf(file = "~/Desktop/Nanopore/Analysis/Differential usage/DEXSeq_DNAJA3.pdf", width = 11.5, height = 8)

plotDEXSeq(object = dxr, geneID = dxr[order(dxr$pvalue), ]$groupID[3], expression = FALSE, splicing = TRUE, displayTranscripts = FALSE, legend = TRUE, color = c("#4040FF", "#FF4162"), cex.axis = 1.2, cex = 2, lwd = 3)

dev.off()

pdf(file = "~/Desktop/Nanopore/Analysis/Differential usage/DEXSeq_TP53TG5.pdf", width = 11.5, height = 8)

plotDEXSeq(object = dxr, geneID = dxr[order(dxr$pvalue), ]$groupID[9], expression = FALSE, splicing = TRUE, displayTranscripts = FALSE, legend = TRUE, color = c("#4040FF", "#FF4162"), cex.axis = 1.2, cex = 2, lwd = 3)

dev.off()

theme_new <- theme_classic(base_size = 24) +
  theme(axis.title.x = element_blank(), 
        legend.position = "bottom", 
        legend.title = element_blank(), 
        text = element_text(family = "Helvetica"))
theme_set(new = theme_new)

samples <- read.csv(file = "~/Desktop/Nanopore/Analysis/Differential expression/colData.csv")
samples$Batch <- factor(x = samples$Batch)
samples$Time <- factor(x = samples$Time)

colnames(x = samples)[c(1, 9)] <- c("sample_id", "condition")

counts <- data.frame(read.table(file = "~/Desktop/Nanopore/Bambu/bambu2/counts_transcript.txt", header = TRUE, sep = "\t"))

colnames(x = counts) <- sub(pattern = "^(.*?)_(.*?)_.*", "\\1_\\2", x = colnames(x = counts))
colnames(x = counts)[1:2] <- c("feature_id", "gene_id")

d <- dmDSdata(counts = counts, samples = samples)
d <- dmFilter(x = d, min_samps_gene_expr = 36, min_samps_feature_expr = 18, min_gene_expr = 10, min_feature_expr = 10)

design <- model.matrix(~ Age + PC1 + PC2 + PC3 + RIN...13:condition + Batch + condition, data = DRIMSeq::samples(x = d))

set.seed(seed = 123)

d <- dmPrecision(x = d, design = design)
d <- dmFit(x = d, design = design)
d <- dmTest(x = d, coef = "conditionTime_2")

res <- DRIMSeq::results(x = d)
res <- res[order(res$pvalue), ]
res$names <- sub(pattern = "\\..*", replacement = "", x = res$gene_id)

Precision   <- plotPrecision(x = d) + labs(x = bquote(Log[10] ~ "of mean expression"), y = bquote(Log[10] ~ "precision")) + theme_classic(base_size = 24) + theme(legend.position = "top", text = element_text(family = "Helvetica"))
PValues     <- plotPValues(x = d, level = "gene") + labs(x = "P-values") + theme_classic(base_size = 24) + theme(text = element_text(family = "Helvetica"))
Proportions1 <- plotProportions(x = d, gene_id = res$gene_id[2], group_variable = "condition", plot_type = "boxplot1") + theme_new + theme(axis.text.x = element_text(size = 14)) + labs(title = "PPM1F")
Proportions2 <- plotProportions(x = d, gene_id = res$gene_id[4], group_variable = "condition", plot_type = "boxplot1") + theme_new + theme(axis.text.x = element_text(size = 14)) + labs(title = "TLK2")
Proportions3 <- plotProportions(x = d, gene_id = res$gene_id[5], group_variable = "condition", plot_type = "boxplot1") + theme_new + theme(axis.text.x = element_text(size = 14)) + labs(title = "MAPK14")

pScreen <- DRIMSeq::results(x = d)$pvalue

names(x = pScreen) <- DRIMSeq::results(x = d)$gene_id

pScreen <- na.omit(object = pScreen)

pConfirmation <- matrix(data = DRIMSeq::results(x = d, level = "feature")$pvalue, ncol = 1)

rownames(x = pConfirmation) <- DRIMSeq::results(x = d, level = "feature")$feature_id

pConfirmation <- na.omit(object = pConfirmation)

tx2gene <- DRIMSeq::results(x = d, level = "feature")[, c("feature_id", "gene_id")]

stageRObj <- stageRTx(pScreen = pScreen, pConfirmation = pConfirmation, pScreenAdjusted = FALSE, tx2gene = tx2gene)
stageRObj <- stageWiseAdjustment(object = stageRObj, method = "dtu", alpha = 0.1, allowNA = FALSE)

SignificantGenes <- getSignificantGenes(object = stageRObj)
SignificantTX    <- getSignificantTx(object = stageRObj)

padj <- getAdjustedPValues(object = stageRObj, order = TRUE, onlySignificantGenes = FALSE)
padj$gene_id <- sub(pattern = "\\..*", replacement = "", x = padj$geneID)

gostres <- gost(query = unique(x = res[startsWith(x = res$names, prefix = "ENS") & res$pvalue < 0.05, ]$names), organism = "hsapiens", ordered_query = TRUE, significant = TRUE)

ggsave(filename = "~/Desktop/Nanopore/Analysis/Differential usage/Precision.pdf", plot = Precision, width = 11, height = 8.5)
ggsave(filename = "~/Desktop/Nanopore/Analysis/Differential usage/PValues.pdf", plot = PValues, width = 11, height = 8.5)
ggsave(filename = "~/Desktop/Nanopore/Analysis/Differential usage/Proportions1.pdf", plot = Proportions1, width = 11, height = 8.5)
ggsave(filename = "~/Desktop/Nanopore/Analysis/Differential usage/Proportions2.pdf", plot = Proportions2, width = 11, height = 8.5)
ggsave(filename = "~/Desktop/Nanopore/Analysis/Differential usage/Proportions3.pdf", plot = Proportions3, width = 11, height = 8.5)

write.csv(x = as.data.frame(x = gostres$result)[, -14], file = "~/Desktop/Nanopore/Analysis/Differential usage/gostres_DRIMSeq.csv", row.names = FALSE)
write.csv(x = res, file = "~/Desktop/Nanopore/Analysis/Differential usage/DRIMSeq.csv", row.names = FALSE)