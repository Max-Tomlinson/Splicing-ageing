library(data.table)
library(ggplot2)
library(ggsci)
library(gprofiler2)

res <- NULL

counts_transcript    <- read.table(file = "~/Desktop/Nanopore/Bambu/bambu2/counts_transcript.txt", header = TRUE, sep = '\t')

custom_bg_transcript <- unique(x = sub(pattern = "\\..*", replacement = "", x = counts_transcript$TXNAME))
custom_bg_gene       <- unique(x = sub(pattern = "\\..*", replacement = "", x = counts_transcript$GENEID))

p_values <- list()

for (i in list.files(path = "~/Desktop/Nanopore/Analysis/glm", pattern = "results.txt", full.names = TRUE)) {
  
  w1 <- fread(input = i, header = TRUE, sep = "\t", fill = TRUE)
  w1$transcript <- sub(pattern = ".*;", replacement = "", x = w1$genes)
  w1$transcript <- sub(pattern = "\\..*", replacement = "", x = w1$transcript)
  w1$genes <- sub(pattern = ";.*", replacement = "", x = w1$genes)
  w1$genes <- sub(pattern = "\\..*", replacement = "", x = w1$genes)
  
  w2 <- w1[which(x = w1$adjusted_p_values < 0.05), ]
  w2 <- w2[order(w2$raw_p_values), ]
  
  res <- rbind(res, cbind(Events = nrow(x = w1), `FDR5%` = nrow(x = w2), Up = nrow(x = w2[which(x = w2$mean_difference > 0), ]), Down = nrow(x = w2[which(x = w2$mean_difference < 0), ]), Genes = length(x = unique(x = w2$genes))))
  
  p_values[[basename(path = i)]] <- w1$raw_p_values
  
  if (!grepl(pattern = "isoform|gini|laplace", x = i)){

    gostres <- gost(query = unique(x = w2$genes), organism = "hsapiens", ordered_query = TRUE, significant = TRUE, correction_method = "gSCS", domain_scope = "custom", custom_bg = custom_bg_gene)

    write.csv(x = as.data.frame(x = gostres$result)[, -14], file = paste0("~/Desktop/Nanopore/Analysis/glm/gostres/", gsub(pattern = "fitted_results.txt", replacement = "gostres.csv", x = basename(path = i))), row.names = FALSE)

  } else if (grepl(pattern = "isoform", x = i)){

    gostres <- gost(query = unique(x = w2$transcript), organism = "hsapiens", ordered_query = TRUE, significant = TRUE, correction_method = "gSCS", domain_scope = "custom", custom_bg = custom_bg_transcript)

    write.csv(x = as.data.frame(x = gostres$result)[, -14], file = paste0("~/Desktop/Nanopore/Analysis/glm/gostres/", gsub(pattern = "fitted_results.txt", replacement = "gostres.csv", x = basename(path = i))), row.names = FALSE)

  } else {

    up   <- w2[which(x = w2$mean_difference > 0), ]
    down <- w2[which(x = w2$mean_difference < 0), ]

    gostres_up   <- gost(query = unique(x = up$genes), organism = "hsapiens", ordered_query = TRUE, significant = TRUE, correction_method = "gSCS", domain_scope = "custom", custom_bg = custom_bg_gene)
    gostres_down <- gost(query = unique(x = down$genes), organism = "hsapiens", ordered_query = TRUE, significant = TRUE, correction_method = "gSCS", domain_scope = "custom", custom_bg = custom_bg_gene)

    up_terms   <- gostres_up$result$term_name
    down_terms <- gostres_down$result$term_name

    unique_up_terms   <- setdiff(up_terms, down_terms)
    unique_down_terms <- setdiff(down_terms, up_terms)

    gostres_up    <- gostres_up$result[gostres_up$result$term_name %in% unique_up_terms, ]
    gostres_down  <- gostres_down$result[gostres_down$result$term_name %in% unique_down_terms, ]

    write.csv(x = gostres_up[, -14], file = paste0("~/Desktop/Nanopore/Analysis/glm/gostres/", gsub(pattern = "fitted_results.txt", replacement = "gostres_up.csv", x = basename(path = i))), row.names = FALSE)
    write.csv(x = gostres_down[, -14], file = paste0("~/Desktop/Nanopore/Analysis/glm/gostres/", gsub(pattern = "fitted_results.txt", replacement = "gostres_down.csv", x = basename(path = i))), row.names = FALSE)

  }
}

res <- as.data.frame(x = res)

names <- gsub(pattern = "bambu_", replacement = "", x = list.files(path = "~/Desktop/Nanopore/Analysis/glm", pattern = "results.txt"))
names <- gsub(pattern = "_fitted_results.txt", replacement = "", x = names)

rownames(x = res) <- names

write.csv(x = res, file = "~/Desktop/Nanopore/Analysis/glm/results.csv", quote = FALSE, row.names = TRUE)

data <- data.table(file = character(), p_value = numeric())

for (i in seq_along(along.with = p_values)) {
  file_name <- basename(path = names(x = p_values[i]))
  data <- rbindlist(l = list(data, data.table(file = file_name, p_value = p_values[[i]])), use.names = TRUE)
}

data$file <- gsub(pattern = "bambu_", replacement = "", x = data$file)
data$file <- gsub(pattern = "_fitted_results.txt", replacement = "", x = data$file)
data <- data[!which(x = data$file == "isoform"), ]

plot1 <- ggplot(data = data[!data$file %in% c("gini_index", "laplace_entropy"), ], mapping = aes(x = p_value, fill = factor(x = file))) +
  geom_histogram(alpha = 0.5, position = "identity") +
  labs(x = "P-value", y = "Frequency") +
  scale_fill_futurama() +
  theme_classic(base_size = 24) +
  theme(legend.position = "bottom", legend.title = element_blank(), text = element_text(family = "Helvetica"))

plot2 <- ggplot(data = data[data$file %in% c("gini_index", "laplace_entropy"), ], mapping = aes(x = p_value, fill = factor(x = file))) +
  geom_histogram(alpha = 0.5, position = "identity") +
  labs(x = "P-value", y = "Frequency") +
  scale_fill_futurama(labels = c("Gini index", "Laplace entropy")) +
  theme_classic(base_size = 24) +
  theme(legend.position = "bottom", legend.title = element_blank(), text = element_text(family = "Helvetica"))

ggsave(filename = "~/Desktop/Nanopore/Analysis/glm/pvalues1.pdf", plot = plot1, width = 11, height = 8.5, useDingbats = FALSE)
ggsave(filename = "~/Desktop/Nanopore/Analysis/glm/pvalues2.pdf", plot = plot2, width = 11, height = 8.5, useDingbats = FALSE)