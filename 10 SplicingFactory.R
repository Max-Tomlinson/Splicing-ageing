library(ggplot2)
library(gprofiler2)
library(ggsci)
library(SplicingFactory)
library(SummarizedExperiment)
library(tidyr)

theme_new <- theme_classic(base_size = 24) +
  theme(legend.position = "bottom", 
        legend.title = element_blank(), 
        text = element_text(family = "Helvetica"))
theme_set(new = theme_new)

counts_transcript <- read.table(file = "~/Desktop/Nanopore/Bambu/bambu2/counts_transcript.txt", header = TRUE, sep = '\t')

genes <- counts_transcript[, 2]

readcounts <- counts_transcript[, -c(1:2)]

laplace_entropy <- calculate_diversity(x = readcounts, genes = genes, method = "laplace", norm = TRUE, verbose = TRUE)

gini_index <- calculate_diversity(x = readcounts, genes = genes, method = "gini", verbose = TRUE)

laplace_data <- cbind(assay(x = laplace_entropy), Gene = rowData(x = laplace_entropy)$genes)
laplace_data <- pivot_longer(data = laplace_data, cols = -Gene, names_to = "sample", values_to = "entropy")
laplace_data$sample_type <- apply(X = laplace_data[, 2], MARGIN = 1, FUN = function(x) ifelse(test = grepl(pattern = "^sample_0[0-9]|^sample_1[0-9]|^sample_2[0-5]", x), yes = "Time 1", no = "Time 2"))
laplace_data <- drop_na(data = laplace_data)
laplace_data$Gene <- paste0(laplace_data$Gene, "_", laplace_data$sample_type)
laplace_data$diversity <-  "Laplace entropy"

gini_data <- cbind(assay(x = gini_index), Gene = rowData(x = gini_index)$genes)
gini_data <- pivot_longer(data = gini_data, cols = -Gene, names_to = "sample", values_to = "gini")
gini_data$sample_type <- apply(X = gini_data[, 2], MARGIN = 1, FUN = function(x) ifelse(test = grepl(pattern = "^sample_0[0-9]|^sample_1[0-9]|^sample_2[0-5]", x), yes = "Time 1", no = "Time 2"))
gini_data <- drop_na(data = gini_data)
gini_data$Gene <- paste0(gini_data$Gene, "_", gini_data$sample_type)
gini_data$diversity <-  "Gini index"

laplace_entropy_mean <- aggregate(x = laplace_data$entropy, by = list(laplace_data$Gene), mean)
colnames(x = laplace_entropy_mean)[2] <- "mean_entropy"
laplace_entropy_mean <- as_tibble(x = laplace_entropy_mean)
laplace_entropy_mean$sample_type <- apply(X = laplace_entropy_mean[, 1], MARGIN = 1, FUN = function(x) ifelse(test = grepl(pattern = "Time 1", x), yes = "Time 1", no = "Time 2"))
laplace_entropy_mean$diversity <-  "Laplace entropy"

gini_mean <- aggregate(x = gini_data$gini, by = list(gini_data$Gene), mean)
colnames(x = gini_mean)[2] <- "mean_gini"
gini_mean <- as_tibble(x = gini_mean)
gini_mean$sample_type <- apply(X = gini_mean[, 1], MARGIN = 1, FUN = function(x) ifelse(test = grepl(pattern = "Time 1", x), yes = "Time 1", no = "Time 2"))
gini_mean$diversity <-  "Gini index"

colData(x = laplace_entropy) <- cbind(colData(x = laplace_entropy), sample_type = ifelse(test = grepl(pattern = "^sample_0[0-9]|^sample_1[0-9]|^sample_2[0-5]", x = laplace_entropy$samples), yes = "Time 1", no = "Time 2"))
colData(x = gini_index) <- cbind(colData(x = gini_index), sample_type = ifelse(test = grepl(pattern = "^sample_0[0-9]|^sample_1[0-9]|^sample_2[0-5]", x = gini_index$samples), yes = "Time 1", no = "Time 2"))

entropy_significance <- calculate_difference(x = laplace_entropy, samples = "sample_type", control = "Time 1", method = "mean", test = "wilcoxon", verbose = TRUE)
entropy_significance$label <- apply(X = entropy_significance[, c(4, 6)], MARGIN = 1, FUN = function(x) ifelse(test = x[1] >= 0.1 & x[2] < 0.05, yes = "Increased", no = ifelse(test = x[1] <= -0.1 & x[2] < 0.05, yes = "Reduced", no = "Non-significant")))
entropy_significance$mean <- apply(X = entropy_significance[, c(2, 3)], MARGIN = 1, FUN = function(x) (x[1] + x[2]) / 2)
entropy_significance$genes <- sub(pattern = "\\..*", replacement = "", x = entropy_significance$genes)
entropy_significance <- entropy_significance[order(entropy_significance$raw_p_values), ]

sig_genes <- entropy_significance[which(x = entropy_significance$raw_p_values < 0.05), ]
sig_genes_up <- sig_genes[which(x = sig_genes$mean_difference > 0), ]
sig_genes_down <- sig_genes[which(x = sig_genes$mean_difference < 0), ]

gostres <- gost(query = sig_genes$genes, organism = "hsapiens", ordered_query = TRUE, significant = TRUE, correction_method = "gSCS", domain_scope = "custom", custom_bg = entropy_significance$genes)
gostres_up <- gost(query = sig_genes_up$genes, organism = "hsapiens", ordered_query = TRUE, significant = TRUE, correction_method = "gSCS", domain_scope = "custom", custom_bg = entropy_significance$genes)
gostres_down <- gost(query = sig_genes_down$genes, organism = "hsapiens", ordered_query = TRUE, significant = TRUE, correction_method = "gSCS", domain_scope = "custom", custom_bg = entropy_significance$genes)

write.csv(x = sig_genes, "~/Desktop/Nanopore/Analysis/Transcript diversity/sig_genes_laplace.csv", row.names = FALSE)
write.csv(x = sig_genes_up, "~/Desktop/Nanopore/Analysis/Transcript diversity/sig_genes_up_laplace.csv", row.names = FALSE)
write.csv(x = sig_genes_down, "~/Desktop/Nanopore/Analysis/Transcript diversity/sig_genes_down_laplace.csv", row.names = FALSE)

write.csv(x = as.data.frame(gostres$result)[, -14], "~/Desktop/Nanopore/Analysis/Transcript diversity/gostres_laplace.csv", row.names = FALSE)
write.csv(x = as.data.frame(gostres_up$result)[, -14], "~/Desktop/Nanopore/Analysis/Transcript diversity/gostres_up_laplace.csv", row.names = FALSE)
write.csv(x = as.data.frame(gostres_down$result)[, -14], "~/Desktop/Nanopore/Analysis/Transcript diversity/gostres_down_laplace.csv", row.names = FALSE)

gini_significance <- calculate_difference(x = gini_index, samples = "sample_type", control = "Time 1", method = "mean", test = "wilcoxon", verbose = TRUE)
gini_significance$label <- apply(X = gini_significance[, c(4, 6)], MARGIN = 1, FUN = function(x) ifelse(test = x[1] >= 0.1 & x[2] < 0.05, yes = "Increased", no = ifelse(test = x[1] <= -0.1 & x[2] < 0.05, yes = "Reduced", no = "Non-significant")))
gini_significance$mean <- apply(X = gini_significance[, c(2, 3)], MARGIN = 1, FUN = function(x) (x[1] + x[2]) / 2)
gini_significance$genes <- sub(pattern = "\\..*", replacement = "", x = gini_significance$genes)
gini_significance <- gini_significance[order(gini_significance$raw_p_values), ]

sig_genes <- gini_significance[which(x = gini_significance$raw_p_values < 0.05), ]
sig_genes_up <- sig_genes[which(x = sig_genes$mean_difference > 0), ]
sig_genes_down <- sig_genes[which(x = sig_genes$mean_difference < 0), ]

gostres <- gost(query = sig_genes$genes, organism = "hsapiens", ordered_query = TRUE, significant = TRUE, correction_method = "gSCS", domain_scope = "custom", custom_bg = entropy_significance$genes)
gostres_up <- gost(query = sig_genes_up$genes, organism = "hsapiens", ordered_query = TRUE, significant = TRUE, correction_method = "gSCS", domain_scope = "custom", custom_bg = entropy_significance$genes)
gostres_down <- gost(query = sig_genes_down$genes, organism = "hsapiens", ordered_query = TRUE, significant = TRUE, correction_method = "gSCS", domain_scope = "custom", custom_bg = entropy_significance$genes)

write.csv(x = sig_genes, "~/Desktop/Nanopore/Analysis/Transcript diversity/sig_genes_gini.csv", row.names = FALSE)
write.csv(x = sig_genes_up, "~/Desktop/Nanopore/Analysis/Transcript diversity/sig_genes_up_gini.csv", row.names = FALSE)
write.csv(x = sig_genes_down, "~/Desktop/Nanopore/Analysis/Transcript diversity/sig_genes_down_gini.csv", row.names = FALSE)

write.csv(x = as.data.frame(gostres$result)[, -14], "~/Desktop/Nanopore/Analysis/Transcript diversity/gostres_gini.csv", row.names = FALSE)
write.csv(x = as.data.frame(gostres_up$result)[, -14], "~/Desktop/Nanopore/Analysis/Transcript diversity/gostres_up_gini.csv", row.names = FALSE)
write.csv(x = as.data.frame(gostres_down$result)[, -14], "~/Desktop/Nanopore/Analysis/Transcript diversity/gostres_down_gini.csv", row.names = FALSE)

ggplot() +
  facet_grid(. ~ sample_type) +
  geom_density(data = laplace_data, alpha = 0.3, mapping = aes(x = entropy, group = sample, colour = diversity)) +
  geom_density(data = gini_data, alpha = 0.3, mapping = aes(x = gini, group = sample, colour = diversity)) +
  guides(colour = "none") +
  labs(x = "Diversity values", y = "Density") +
  scale_colour_futurama()

ggsave(filename = "~/Desktop/Nanopore/Analysis/Transcript diversity/Density_diversity.pdf", last_plot(), width = 11, height = 8.5)

ggplot() +
  geom_density(data = laplace_data, mapping = aes(x = entropy, group = sample_type, colour = sample_type)) +
  labs(x = "Laplace entropy", y = "Density") +
  scale_color_manual(name = "Time", breaks = c("Time 1", "Time 2"), labels = c("1", "2"), values = c("Time 1" = "#475da0", "Time 2" = "#ba1f21")) +
  theme(legend.position = "top")

ggsave(filename = "~/Desktop/Nanopore/Analysis/Transcript diversity/Density_laplace.pdf", last_plot(), width = 11, height = 8.5)

ggplot() +
  geom_density(data = gini_data, mapping = aes(x = gini, group = sample_type, colour = sample_type)) +
  labs(x = "Gini index", y = "Density") +
  scale_colour_futurama(name = "Time", breaks = c("Time 1", "Time 2"), labels = c("1", "2"))

ggsave(filename = "~/Desktop/Nanopore/Analysis/Transcript diversity/Density_gini.pdf", last_plot(), width = 11, height = 8.5)

ggplot() +
  coord_flip() +
  geom_violin(data = laplace_entropy_mean, aes(x = sample_type, y = mean_entropy, fill = diversity), alpha = 0.6) +
  geom_violin(data = gini_mean, aes(x = sample_type, y = mean_gini, fill = diversity), alpha = 0.6) +
  labs(x = NULL, y = "Diversity") +
  scale_fill_futurama(name = "Diversity")

ggsave(filename = "~/Desktop/Nanopore/Analysis/Transcript diversity/Violin_diversity.pdf", last_plot(), width = 11, height = 8.5)

entropy_significance <- entropy_significance[!is.na(x = entropy_significance$raw_p_values), ]

ggplot(data = entropy_significance, mapping = aes(x = mean, y = mean_difference, colour = label)) +
  geom_point(colour = "grey", size = 0.5) +
  geom_point(data = entropy_significance[!entropy_significance$label == "Non-significant", ], size = 1) +
  labs(x = "Mean Laplace entropy", y = "Mean difference") + 
  scale_color_manual(values = c("Increased" = "#ba1f21", "Reduced" = "#475da0"), name = NULL)

ggsave(filename = "~/Desktop/Nanopore/Analysis/Transcript diversity/MA_entropy.pdf", last_plot(), width = 11, height = 8.5)

ggplot(data = entropy_significance, mapping = aes(x = mean_difference, y = -log10(x = raw_p_values), colour = label)) +
  geom_point(colour = "grey", size = 0.5) +
  geom_point(data = entropy_significance[!entropy_significance$label == "Non-significant", ], size = 1) +
  geom_hline(yintercept = -log10(x = 0.05), color = "red") +
  geom_vline(xintercept = c(0.1, -0.1)) +
  labs(x = "Mean difference in Laplace entropy between time points", y = expression(-Log[10]~italic(P))) +
  scale_color_manual(values = c("Increased" = "#ba1f21", "Reduced" = "#475da0"), guide = "none")

ggsave(filename = "~/Desktop/Nanopore/Analysis/Transcript diversity/MA_entropy2.pdf", last_plot(), width = 11, height = 8.5)

gini_significance <- gini_significance[!is.na(x = gini_significance$raw_p_values), ]

ggplot(data = gini_significance, mapping = aes(x = mean, y = mean_difference, colour = label)) +
  geom_point(colour = "grey", size = 0.5) +
  geom_point(data = gini_significance[!gini_significance$label == "Non-significant", ], size = 1) +
  labs(x = "Mean Gini index", y = "Mean difference") + 
  scale_color_manual(values = c("Increased" = "#ba1f21", "Reduced" = "#475da0"), name = NULL)

ggsave(filename = "~/Desktop/Nanopore/Analysis/Transcript diversity/MA_gini.pdf", last_plot(), width = 11, height = 8.5)

ggplot(data = gini_significance, mapping = aes(x = mean_difference, y = -log10(x = raw_p_values), colour = label)) +
  geom_point(colour = "grey", size = 0.5) +
  geom_point(data = gini_significance[!gini_significance$label == "Non-significant", ], size = 1) +
  geom_hline(yintercept = -log10(x = 0.05), color = "red") +
  geom_vline(xintercept = c(0.1, -0.1)) +
  labs(x = "Mean difference in Gini index between time points", y = expression(-Log[10]~italic(P))) +
  scale_color_manual(values = c("Increased" = "#ba1f21", "Reduced" = "#475da0"), guide = "none")

ggsave(filename = "~/Desktop/Nanopore/Analysis/Transcript diversity/MA_gini2.pdf", last_plot(), width = 11, height = 8.5)
