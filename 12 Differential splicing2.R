library(ggplot2)
library(gprofiler2)
library(ggsci)
library(tidyr)

theme_new <- theme_classic(base_size = 24) +
  theme(legend.position = "none", 
        legend.title = element_blank(), 
        text = element_text(family = "Helvetica"))
theme_set(new = theme_new)

laplace_data <- read.table(file = "~/Desktop/Nanopore/Analysis/glm/laplace_entropy_fitted.txt", header = TRUE, sep = '\t', row.names = 1)
laplace_data <- cbind(laplace_data, Gene = rownames(x = laplace_data))
laplace_data <- pivot_longer(data = laplace_data, cols = -Gene, names_to = "sample", values_to = "entropy")
laplace_data$sample_type <- apply(X = laplace_data[, 2], MARGIN = 1, FUN = function(x) ifelse(test = grepl(pattern = "^sample_0[0-9]|^sample_1[0-9]|^sample_2[0-5]", x), yes = "Time 1", no = "Time 2"))
laplace_data <- drop_na(data = laplace_data)
laplace_data$Gene <- paste0(laplace_data$Gene, "_", laplace_data$sample_type)
laplace_data$diversity <-  "Laplace entropy"

gini_data <- read.table(file = "~/Desktop/Nanopore/Analysis/glm/gini_index_fitted.txt", header = TRUE, sep = '\t', row.names = 1)
gini_data <- cbind(gini_data, Gene = rownames(x = gini_data))
gini_data <- pivot_longer(data = gini_data, cols = -Gene, names_to = "sample", values_to = "gini")
gini_data$sample_type <- apply(X = gini_data[, 2], MARGIN = 1, FUN = function(x) ifelse(test = grepl(pattern = "^sample_0[0-9]|^sample_1[0-9]|^sample_2[0-5]", x), yes = "Time 1", no = "Time 2"))
gini_data <- drop_na(data = gini_data)
gini_data$Gene <- paste0(gini_data$Gene, "_", gini_data$sample_type)
gini_data$diversity <-  "Gini index"

laplace_entropy_mean <- aggregate(x = laplace_data$entropy, by = list(laplace_data$Gene), mean)
colnames(x = laplace_entropy_mean)[2] <- "mean_entropy"
laplace_entropy_mean <- as_tibble(x = laplace_entropy_mean)
laplace_entropy_mean$sample_type <- apply(X = laplace_entropy_mean[, 1], MARGIN = 1, FUN = function(x) ifelse(test = grepl(pattern = "Time 1", x), yes = "Time 1", no = "Time 2"))

gini_mean <- aggregate(x = gini_data$gini, by = list(gini_data$Gene), mean)
colnames(x = gini_mean)[2] <- "mean_gini"
gini_mean <- as_tibble(x = gini_mean)
gini_mean$sample_type <- apply(X = gini_mean[, 1], MARGIN = 1, FUN = function(x) ifelse(test = grepl(pattern = "Time 1", x), yes = "Time 1", no = "Time 2"))

entropy_significance <- read.table(file = "~/Desktop/Nanopore/Analysis/glm/laplace_entropy_fitted_results.txt", header = TRUE, sep = '\t', row.names = 1)
entropy_significance$label <- apply(X = entropy_significance[, c(3, 6)], MARGIN = 1, FUN = function(x) ifelse(test = x[1] >= 0 & x[2] < 0.05, yes = "Increased Laplace entropy", no = ifelse(test = x[1] <= -0 & x[2] < 0.05, yes = "Reduced Laplace entropy", no = "Non-significant")))
entropy_significance$mean <- apply(X = entropy_significance[, c(1, 2)], MARGIN = 1, FUN = function(x) (x[1] + x[2]) / 2)
entropy_significance$genes <- sub(pattern = "\\..*", replacement = "", x = rownames(x = entropy_significance))
entropy_significance <- entropy_significance[order(entropy_significance$raw_p_values), ]

sig_genes <- entropy_significance[which(x = entropy_significance$adjusted_p_values < 0.05), ]

write.csv(x = sig_genes[which(x = sig_genes$mean_difference > 0), ], "~/Desktop/Nanopore/Analysis/glm/sig_genes_up_laplace.csv", row.names = FALSE)
write.csv(x = sig_genes[which(x = sig_genes$mean_difference < 0), ], "~/Desktop/Nanopore/Analysis/glm/sig_genes_down_laplace.csv", row.names = FALSE)

gini_significance <- read.table(file = "~/Desktop/Nanopore/Analysis/glm/gini_index_fitted_results.txt", header = TRUE, sep = '\t', row.names = 1)
gini_significance$label <- apply(X = gini_significance[, c(3, 6)], MARGIN = 1, FUN = function(x) ifelse(test = x[1] >= 0 & x[2] < 0.05, yes = "Increased Gini index", no = ifelse(test = x[1] <= -0 & x[2] < 0.05, yes = "Reduced Gini index", no = "Non-significant")))
gini_significance$mean <- apply(X = gini_significance[, c(1, 2)], MARGIN = 1, FUN = function(x) (x[1] + x[2]) / 2)
gini_significance$genes <- sub(pattern = "\\..*", replacement = "", x = rownames(x = gini_significance))
gini_significance <- gini_significance[order(gini_significance$raw_p_values), ]

ggplot(data = laplace_entropy_mean, mapping = aes(x = sample_type, y = mean_entropy, fill = sample_type)) +
  coord_flip() +
  geom_violin() +
  geom_boxplot(width = 0.25, fill = "white") +
  guides(fill = "none") +
  labs(x = NULL, y = "Diversity") +
  scale_fill_futurama()

ggsave(filename = "~/Desktop/Nanopore/Analysis/glm/Violin_diversity_fitted.pdf", last_plot(), width = 11, height = 8.5)

ggplot(data = entropy_significance, mapping = aes(x = mean, y = mean_difference, colour = label)) +
  geom_point(colour = "grey", size = 0.5) +
  geom_point(data = entropy_significance[!entropy_significance$label == "Non-significant", ], size = 1) +
  labs(x = "Mean Laplace entropy", y = "Mean difference") + 
  scale_color_manual(values = c("Increased Laplace entropy" = "#ba1f21", "Reduced Laplace entropy" = "#475da0"), name = NULL)

ggsave(filename = "~/Desktop/Nanopore/Analysis/glm/MA_entropy_fitted.pdf", last_plot(), width = 11, height = 8.5)

ggplot(data = gini_significance, mapping = aes(x = mean, y = mean_difference, colour = label)) +
  geom_point(colour = "grey", size = 0.5) +
  geom_point(data = gini_significance[!gini_significance$label == "Non-significant", ], size = 1) +
  labs(x = "Mean Gini index", y = "Mean difference") + 
  scale_color_manual(values = c("Increased Gini index" = "#ba1f21", "Reduced Gini index" = "#475da0"), name = NULL)

ggsave(filename = "~/Desktop/Nanopore/Analysis/glm/MA_gini_fitted.pdf", last_plot(), width = 11, height = 8.5)
