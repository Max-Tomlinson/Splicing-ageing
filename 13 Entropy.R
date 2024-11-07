library(entropy)
library(readxl)

cts <- read.table(file = "~/Desktop/Nanopore/Bambu/bambu2/counts_gene.txt", header = TRUE, sep = '\t')
cts$GENEID <- sub(pattern = "\\..*$", replacement = "", x = cts$GENEID)

genes_up   <- read.csv(file = "~/Desktop/Nanopore/Analysis/glm/sig_genes_up_laplace.csv", stringsAsFactors = FALSE)[1:1000, ]
genes_down <- read.csv(file = "~/Desktop/Nanopore/Analysis/glm/sig_genes_down_laplace.csv", stringsAsFactors = FALSE)[1:1000, ]

cts_up   <- cts[cts$GENEID %in% genes_up$genes, ]
cts_down <- cts[cts$GENEID %in% genes_down$genes, ]

entropy_up <- lapply(X = cts_up[, -1], FUN = function(column) entropy(y = column, method = "Laplace"))
entropy_up <- as.data.frame(do.call(what = cbind, args = entropy_up), row.names = "Laplace")
entropy_down <- lapply(X = cts_down[, -1], FUN = function(column) entropy(y = column, method = "Laplace"))
entropy_down <- as.data.frame(do.call(what = cbind, args = entropy_down), row.names = "Laplace")
entropy <- merge(x = t(x = entropy_up), y = t(x = entropy_down), by = 0, suffixes = c("_up", "_down"))
entropy$Laplace_diff <- entropy$Laplace_up - entropy$Laplace_down

covar <- read.table(file = "~/Desktop/Nanopore/Analysis/Differential expression/dCell.csv", header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE)

extraction_data <- read_xlsx(path = "~/Documents/KCL/PhD/TwinsUK/Data/Nanopore/Extraction data/RNA_extraction_data_updated1[1].xlsx", skip = 1)
extraction_data <- extraction_data[-c(12, 14, 32, 34), c("RIN...13")]

data <- cbind.data.frame(covar, RIN = extraction_data$RIN...13, entropy[, -1])

sample_ids <- colnames(x = cts[, -1])

rownames(x = data) <- sample_ids

fitted <- NULL

for (i in 9:ncol(x = data)){
  
  model <- lm(formula = as.formula(object = paste(colnames(x = data)[i], "~ Age + PC1 + PC2 + PC3 + RIN:factor(x = Sample)")), data = data)

  fitted <- rbind(fitted, fitted(object = model))
}

ma_ids <- sample_ids[c(1:8, 19:26)]
oa_ids <- sample_ids[c(9:18, 27:36)]
tp1_ids <- sample_ids[c(1:18)]
tp2_ids <- sample_ids[c(19:36)]

ma_tp1_ids <- intersect(x = ma_ids, y = tp1_ids)
ma_tp2_ids <- intersect(x = ma_ids, y = tp2_ids)
oa_tp1_ids <- intersect(x = oa_ids, y = tp1_ids)
oa_tp2_ids <- intersect(x = oa_ids, y = tp2_ids)

perform_test <- function(df, ids1, ids2, paired = TRUE, alternative = "two.sided") {
  test_results <- apply(X = df, MARGIN = 1, FUN = function(row) {
    group1_numeric <- as.numeric(x = unlist(x = row[ids1]))
    group2_numeric <- as.numeric(x = unlist(x = row[ids2]))
    test_result <- t.test(x = group1_numeric, y = group2_numeric, paired = paired, alternative = alternative)
    test_result <- list(p.value = test_result$p.value, mean_difference = mean(x = group1_numeric) - mean(x = group2_numeric))
    return(test_result)
  })
  
  test_summaries <- lapply(X = test_results, FUN = function(t) { data.frame(mean_difference = signif(x = t$mean_difference, digits = 3), p.value = signif(x = t$p.value, digits = 3))
  })
  
  test_summaries_df <- do.call(what = rbind, args = test_summaries)
  return(test_summaries_df)
}

tp1_vs_tp2 <- perform_test(df = fitted, ids1 = tp2_ids, ids2 = tp1_ids)
tp1_vs_tp2_ma <- perform_test(df = fitted, ids1 = ma_tp2_ids, ids2 = ma_tp1_ids)
tp1_vs_tp2_oa <- perform_test(df = fitted, ids1 = oa_tp2_ids, ids2 = oa_tp1_ids)

results <- cbind.data.frame(tp1_vs_tp2, tp1_vs_tp2_ma, tp1_vs_tp2_oa)

rownames(x = results) <- colnames(x = data)[9:11]
colnames(x = results)[seq(from = 1, to = 6, by = 2)] <- c("T2_v_T1", "T2_v_T1_MA", "T2_v_T1_OA")

write.csv(x = results, file = "~/Desktop/Nanopore/Analysis/Transcript diversity/entropy.csv", row.names = TRUE)