library(DESeq2)
library(EnhancedVolcano)
library(gprofiler2)
library(IHW)
library(qvalue)
library(readxl)

colData <- read.csv(file = "~/Desktop/Nanopore/Analysis/Differential expression/dCell.csv", row.names = 1)

rownames(x = colData) <- sub(pattern = "_full_length_output.sorted.bam", replacement = "", rownames(x = colData))

extraction_data <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Nanopore/Extraction data/RNA_extraction_data_updated1[1].xlsx", skip = 1)

colData <- cbind(colData, extraction_data[-c(12, 14, 32, 34), c("RIN...13")])
colData <- cbind.data.frame(scale(x = colData[, c(3:6, 8)]), Batch = colData$Batch)
colData$Batch <- factor(x = colData$Batch)
colData$Group <- factor(x = c(rep(x = "Middle_aged", 8), rep(x = "Older_aged", 10), rep(x = "Middle_aged", 8), rep(x = "Older_aged", 10)))
colData$Time <- factor(x = c(rep(x = "Time_1", 18), rep(x = "Time_2", 18)))

countData <- read.csv(file = "~/Desktop/Nanopore/Bambu/bambu2/counts_gene.txt", sep = "\t", row.names = 1)

colnames(x = countData) <- sub(pattern = "^(.*?)_(.*?)_.*", replacement = "\\1_\\2", x = colnames(x = countData))

design <- as.formula("~ Age + PC1 + PC2 + PC3 + RIN...13:Time + Batch + Time")

dds <- DESeqDataSetFromMatrix(countData = round(x = countData), colData = colData, design = design)

keep <- rowSums(x = counts(object = dds) >= 10) >= 18

dds <- dds[keep,]
dds <- DESeq(object = dds)

res1 <- DESeq2::results(object = dds, filterFun = ihw)
res1 <- as.data.frame(x = res1[order(res1$padj), ])
res1$gene_id <- sub(pattern = "\\..*", replacement = "", x = rownames(x = res1))

gene_ids <- gconvert(query = res1$gene_id, organism = "hsapiens", mthreshold = 1, filter_na = TRUE)

res1$names <- ifelse(test = !is.na(x = match(x = res1$gene_id, table = gene_ids$input)), yes = gene_ids$name[match(x = res1$gene_id, table = gene_ids$input)], no = res1$gene_id)
res1$names[c(9, 16)] <- c("TEC", "SNORA25")

dds <- DESeqDataSetFromMatrix(countData = round(x = countData[, c(1:8, 19:26)]), colData = colData[which(x = colData$Group == "Middle_aged"), ], design = design)

keep <- rowSums(x = counts(object = dds) >= 10) >= 8

dds <- dds[keep,]
dds <- DESeq(object = dds)

res1.1 <- DESeq2::results(object = dds, filterFun = ihw)
res1.1 <- as.data.frame(x = res1.1[order(res1.1$padj), ])
res1.1$gene_id <- sub(pattern = "\\..*", replacement = "", x = rownames(x = res1.1))

gene_ids <- gconvert(query = res1.1$gene_id, organism = "hsapiens", mthreshold = 1, filter_na = TRUE)

res1.1$names <- ifelse(test = !is.na(x = match(x = res1.1$gene_id, table = gene_ids$input)), yes = gene_ids$name[match(x = res1.1$gene_id, table = gene_ids$input)], no = res1.1$gene_id)
#res1.1$names[c()] <- c("")

dds <- DESeqDataSetFromMatrix(countData = round(x = countData[, c(9:18, 27:36)]), colData = colData[which(x = colData$Group == "Older_aged"), ], design = design)

keep <- rowSums(x = counts(object = dds) >= 10) >= 10

dds <- dds[keep,]
dds <- DESeq(object = dds)

res1.2 <- DESeq2::results(object = dds, filterFun = ihw)
res1.2 <- as.data.frame(x = res1.2[order(res1.2$padj), ])
res1.2$gene_id <- sub(pattern = "\\..*", replacement = "", x = rownames(x = res1.2))

gene_ids <- gconvert(query = res1.2$gene_id, organism = "hsapiens", mthreshold = 1, filter_na = TRUE)

res1.2$names <- ifelse(test = !is.na(x = match(x = res1.2$gene_id, table = gene_ids$input)), yes = gene_ids$name[match(x = res1.2$gene_id, table = gene_ids$input)], no = res1.1$gene_id)
#res1.1$names[c()] <- c("")

MultiMuTHER <- read.table(file = "~/Desktop/Nanopore/Analysis/Differential expression/MultiMuTHER_SupplementaryTable2.txt", header = TRUE, sep = "\t")
MultiMuTHER <- merge(x = res1, y = MultiMuTHER, by.x = 9, by.y = 2, suffixes = c('_1', '_2'), sort = FALSE)
MultiMuTHER <- MultiMuTHER[MultiMuTHER$Population.level.Change.over.time == 1, ]

qobj <- qvalue(p = MultiMuTHER$pvalue)

countData <- read.csv(file = "~/Desktop/Nanopore/Bambu/bambu2/counts_transcript.txt", sep = "\t", row.names = 1)[, -1]

colnames(x = countData) <- sub(pattern = "^(.*?)_(.*?)_.*", replacement = "\\1_\\2", x = colnames(x = countData))

dds <- DESeqDataSetFromMatrix(countData = round(x = countData), colData = colData, design = design)

keep <- rowSums(x = counts(object = dds) >= 10) >= 18

dds <- dds[keep,]
dds <- DESeq(object = dds)

res2 <- DESeq2::results(object = dds, filterFun = ihw)
res2 <- as.data.frame(x = res2[order(res2$padj), ])
res2$transcript_id <- sub(pattern = "\\..*", replacement = "", x = rownames(x = res2))

transcript_ids <- gconvert(query = res2$transcript_id, organism = "hsapiens", mthreshold = 1, filter_na = TRUE)

res2$names <- ifelse(test = !is.na(x = match(x = res2$transcript_id, table = transcript_ids$input)), yes = transcript_ids$name[match(x = res2$transcript_id, table = transcript_ids$input)], no = res2$transcript_id)
#res2$names[c()] <- c("")

dds <- DESeqDataSetFromMatrix(countData = round(x = countData[, c(1:8, 19:26)]), colData = colData[which(x = colData$Group == "Middle_aged"), ], design = design)

keep <- rowSums(x = counts(object = dds) >= 10) >= 8

dds <- dds[keep,]
dds <- DESeq(object = dds)

res2.1 <- DESeq2::results(object = dds, filterFun = ihw)
res2.1 <- as.data.frame(x = res2.1[order(res2.1$padj), ])
res2.1$transcript_id <- sub(pattern = "\\..*", replacement = "", x = rownames(x = res2.1))

transcript_ids <- gconvert(query = res2.1$transcript_id, organism = "hsapiens", mthreshold = 1, filter_na = TRUE)

res2.1$names <- ifelse(test = !is.na(x = match(x = res2.1$transcript_id, table = transcript_ids$input)), yes = transcript_ids$name[match(x = res2.1$transcript_id, table = transcript_ids$input)], no = res2.1$transcript_id)
#res2.1$names[c()] <- c("")

dds <- DESeqDataSetFromMatrix(countData = round(x = countData[, c(9:18, 27:36)]), colData = colData[which(x = colData$Group == "Older_aged"), ], design = design)

keep <- rowSums(x = counts(object = dds) >= 10) >= 10

dds <- dds[keep,]
dds <- DESeq(object = dds)

res2.2 <- DESeq2::results(object = dds, filterFun = ihw)
res2.2 <- as.data.frame(x = res2.2[order(res2.2$padj), ])
res2.2$transcript_id <- sub(pattern = "\\..*", replacement = "", x = rownames(x = res2.2))

transcript_ids <- gconvert(query = res2.2$transcript_id, organism = "hsapiens", mthreshold = 1, filter_na = TRUE)

res2.2$names <- ifelse(test = !is.na(x = match(x = res2.2$transcript_id, table = transcript_ids$input)), yes = transcript_ids$name[match(x = res2.2$transcript_id, table = transcript_ids$input)], no = res2.2$transcript_id)
#res2.2$names[c()] <- c("")

query <- res2[which(x = startsWith(x = res2$transcript_id, prefix = "ENS") & res2$padj < 0.1), ]$transcript_id

custom_bg <- res2[which(x = startsWith(x = res2$transcript_id, prefix = "ENS")), ]$transcript_id

gostres <- gost(query = query, organism = "hsapiens", ordered_query = TRUE, domain_scope = "custom", custom_bg = custom_bg)
gostres$result$query <- "Functional enrichment of differentially expressed transcripts"

gostplot <- gostplot(gostres = gostres, interactive = FALSE) + theme_classic(base_size = 24) + labs(y = expression(-Log[10]~italic(P))) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14), axis.title.x = element_blank(), legend.position = "none", text = element_text(family = "Helvetica"))

publish_gostplot(p = gostplot, highlight_terms = c("WP:WP411", "REAC:R-HSA-72172", "CORUM:351", "KEGG:04140", "WP:WP231"), filename = "~/Desktop/Nanopore/Analysis/Differential expression/gostplot.pdf", width = 11, height = 8.5)

pdf(file = "~/Desktop/Nanopore/Analysis/Differential expression/MultiMuTHER.pdf", width = 11.5, height = 8)

hist(x = qobj) + labs(title = "P-value density histogram of DEGs in MultiMuTHER", x = "P-value", y = "Density") + theme_classic(base_size = 24) + theme(legend.position = "bottom", legend.title = element_blank(), text = element_text(family = "Helvetica"))

dev.off()

lab_italics <- paste0("italic('", make.unique(names = res1$names), "')")

keyvals <- ifelse(test = res1$log2FoldChange < 0 & res1$padj < 0.1, yes = '#475da0', no = ifelse(test = res1$log2FoldChange > 0 & res1$padj < 0.1, yes = '#ba1f21', no = 'grey'))

names(x = keyvals)[keyvals == '#ba1f21'] <- 'Up-regulated'
names(x = keyvals)[keyvals == 'grey'] <- 'NS'
names(x = keyvals)[keyvals == '#475da0'] <- 'Down-regulated'

plot1 <- EnhancedVolcano(toptable = res1,
                        lab = lab_italics,
                        x = 'log2FoldChange',
                        y = 'padj',
                        xlim = range(res1$log2FoldChange, na.rm = TRUE),
                        ylim = c(0, max(-log10(x = res1$padj), na.rm = TRUE)),
                        selectLab = lab_italics[1:20],
                        axisLabSize = 24,
                        title = "Differential gene expression between time points",
                        subtitle = NULL,
                        caption = NULL,
                        titleLabSize = 24,
                        pCutoff = 0.1,
                        pCutoffCol = 'padj',
                        FCcutoff = 0,
                        pointSize = c(ifelse(test = abs(x = res1$log2FoldChange) > 0 & res1$padj < 0.1, yes = 2, no = 2/8)),
                        boxedLabels = TRUE,
                        parseLabels = TRUE,
                        colCustom = keyvals,
                        legendPosition = "none",
                        legendLabSize = 24,
                        drawConnectors = TRUE,
                        arrowheads = FALSE,
                        max.overlaps = 30,
                        gridlines.major = FALSE,
                        gridlines.minor = FALSE)

lab_italics <- paste0("italic('", make.unique(names = res2$names), "')")

keyvals <- ifelse(test = res2$log2FoldChange < 0 & res2$padj < 0.1, yes = '#475da0', no = ifelse(test = res2$log2FoldChange > 0 & res2$padj < 0.1, yes = '#ba1f21', no = 'grey'))

names(x = keyvals)[keyvals == '#ba1f21'] <- 'Up-regulated'
names(x = keyvals)[keyvals == 'grey'] <- 'NS'
names(x = keyvals)[keyvals == '#475da0'] <- 'Down-regulated'

plot2 <- EnhancedVolcano(toptable = res2,
                        lab = lab_italics,
                        x = 'log2FoldChange',
                        y = 'padj',
                        xlim = range(res2$log2FoldChange, na.rm = TRUE),
                        ylim = c(0, max(-log10(x = res2$padj), na.rm = TRUE)),
                        selectLab = lab_italics[1:20],
                        axisLabSize = 24,
                        title = "Differential transcript expression between time points",
                        subtitle = NULL,
                        caption = NULL,
                        titleLabSize = 24,
                        pCutoff = 0.1,
                        pCutoffCol = 'padj',
                        FCcutoff = 0,
                        pointSize = c(ifelse(test = abs(x = res2$log2FoldChange) > 0 & res2$padj < 0.1, yes = 2, no = 2/8)),
                        boxedLabels = TRUE,
                        parseLabels = TRUE,
                        colCustom = keyvals,
                        legendPosition = "none",
                        legendLabSize = 24,
                        drawConnectors = TRUE,
                        arrowheads = FALSE,
                        max.overlaps = 30,
                        gridlines.major = FALSE,
                        gridlines.minor = FALSE)

write.csv(x = as.data.frame(x = gostres$result)[, -14], file = "~/Desktop/Nanopore/Analysis/Differential expression/gostres.csv", row.names = FALSE)
write.csv(x = res1, file = "~/Desktop/Nanopore/Analysis/Differential expression/DEGs.csv", row.names = FALSE)
write.csv(x = res1.1, file = "~/Desktop/Nanopore/Analysis/Differential expression/DEGs_MA.csv", row.names = FALSE)
write.csv(x = res1.2, file = "~/Desktop/Nanopore/Analysis/Differential expression/DEGs_OA.csv", row.names = FALSE)
write.csv(x = res2, file = "~/Desktop/Nanopore/Analysis/Differential expression/DETs.csv", row.names = FALSE)
write.csv(x = res2.1, file = "~/Desktop/Nanopore/Analysis/Differential expression/DETs_MA.csv", row.names = FALSE)
write.csv(x = res2.2, file = "~/Desktop/Nanopore/Analysis/Differential expression/DETs_OA.csv", row.names = FALSE)

table <- rbind(c(n = nrow(x = res1), p.05 = nrow(x = res1[which(x = res1$pvalue < 0.05), ]), FDR5 = nrow(x = res1[which(x = res1$padj < 0.05), ]), Up = nrow(x = res1[which(x = res1$padj < 0.05 & res1$log2FoldChange > 0), ]), Down = nrow(x = res1[which(x = res1$padj < 0.05 & res1$log2FoldChange < 0), ])),
               c(n = nrow(x = res1.1), p.05 = nrow(x = res1.1[which(x = res1.1$pvalue < 0.05), ]), FDR5 = nrow(x = res1.1[which(x = res1.1$padj < 0.05), ]), Up = nrow(x = res1.1[which(x = res1.1$padj < 0.05 & res1.1$log2FoldChange > 0), ]), Down = nrow(x = res1.1[which(x = res1.1$padj < 0.05 & res1.1$log2FoldChange < 0), ])),
               c(n = nrow(x = res1.2), p.05 = nrow(x = res1.2[which(x = res1.2$pvalue < 0.05), ]), FDR5 = nrow(x = res1.2[which(x = res1.2$padj < 0.05), ]), Up = nrow(x = res1.2[which(x = res1.2$padj < 0.05 & res1.2$log2FoldChange > 0), ]), Down = nrow(x = res1.2[which(x = res1.2$padj < 0.05 & res1.2$log2FoldChange < 0), ])),
               c(n = nrow(x = res2), p.05 = nrow(x = res2[which(x = res2$pvalue < 0.05), ]), FDR5 = nrow(x = res2[which(x = res2$padj < 0.05), ]), Up = nrow(x = res2[which(x = res2$padj < 0.05 & res2$log2FoldChange > 0), ]), Down = nrow(x = res2[which(x = res2$padj < 0.05 & res2$log2FoldChange < 0), ])),
               c(n = nrow(x = res2.1), p.05 = nrow(x = res2.1[which(x = res2.1$pvalue < 0.05), ]), FDR5 = nrow(x = res2.1[which(x = res2.1$padj < 0.05), ]), Up = nrow(x = res2.1[which(x = res2.1$padj < 0.05 & res2.1$log2FoldChange > 0), ]), Down = nrow(x = res2.1[which(x = res2.1$padj < 0.05 & res2.1$log2FoldChange < 0), ])),
               c(n = nrow(x = res2.2), p.05 = nrow(x = res2.2[which(x = res2.2$pvalue < 0.05), ]), FDR5 = nrow(x = res2.2[which(x = res2.2$padj < 0.05), ]), Up = nrow(x = res2.2[which(x = res2.2$padj < 0.05 & res2.2$log2FoldChange > 0), ]), Down = nrow(x = res2.2[which(x = res2.2$padj < 0.05 & res2.2$log2FoldChange < 0), ])))

rownames(x = table) <- c("Gene_T1_T2", "Gene_T1_T2_MA", "Gene_T1_T2_OA", "TX_T1_T2", "TX_T1_T2_MA", "TX_T1_T2_OA")

write.csv(x = table, file = "~/Desktop/Nanopore/Analysis/Differential expression/DE_summary.csv", row.names = TRUE)

ggsave(filename = "~/Desktop/Nanopore/Analysis/Differential expression/Volcano_gene.pdf", plot = plot1, width = 11, height = 8.5, useDingbats = FALSE)
ggsave(filename = "~/Desktop/Nanopore/Analysis/Differential expression/Volcano_transcript.pdf", plot = plot2, width = 11, height = 8.5, useDingbats = FALSE)
