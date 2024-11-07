library(caret)
library(DeconCell)
library(dplyr)
library(lubridate)
library(matrixStats)
library(tidyr)
library(truncnorm)

count.table <- read.csv(file = "~/Desktop/Nanopore/Bambu/bambu2/counts_gene.txt", sep = "\t", row.names = 1)

rownames(x = count.table) <- sub(pattern = "\\..*", replacement = "", x = rownames(x = count.table))

data("dCell.models")

dCell <- dCell.expProcessing(count.table = count.table)
dCell <- dCell.predict(dCell.exp = dCell, dCell.models = dCell.models)
dCell <- as.data.frame(x = dCell$dCell.prediction)
dCell2 <- t(x = dCell)
dCell2 <- dCell2[, 19:36] - dCell2[, 1:18]
dCell2 <- t(x = dCell2)

pca <- prcomp(x = dCell, center = TRUE, scale. = TRUE)
pca <- as.data.frame(x = pca$x)
pca2 <- prcomp(x = dCell2, center = TRUE, scale. = TRUE)
pca2 <- as.data.frame(x = pca2$x)

pax_selected_samples <- read.table(file = "~/Documents/KCL/PhD/TwinsUK/Data/Nanopore/pax_selected_samples.txt", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
pax_selected_samples <- pax_selected_samples[!(pax_selected_samples$Study_No %in% c(371, 431, 1832, 2621, 32331, 371, 431, 32861, 51901)), c("Study_No", "Age1", "Age2", "SampleDate1", "SampleDate2")]
pax_selected_samples <- pax_selected_samples %>% pivot_longer(cols = starts_with("SampleDate"), names_to = "Sample", values_to = "Date")
pax_selected_samples <- pax_selected_samples %>% mutate(Age = ifelse(test = Sample == "SampleDate1", yes = Age1, no = Age2))
pax_selected_samples$SampleDate <- dmy(pax_selected_samples$Date)

sample_info <- read.csv(file = "~/Desktop/Nanopore/sample_info.csv")
sample_info <- sample_info[sample_info$Batch != 0, ]
sample_info <- sample_info %>% mutate(Sample = ifelse(test = grepl(pattern = "0[4-9]|1[0-9]|2[0-5]$", x = Names), yes = "SampleDate1", no = "SampleDate2"))
sample_info <- merge(x = pax_selected_samples[, c("Study_No", "Sample", "Age")], y = sample_info[, c("Study_No", "Batch", "Names", "Sample")], by = c("Study_No", "Sample"))
sample_info <- sample_info %>% arrange(Names)

write.csv(x = cbind(sample_info[, 1:3], pca[, 1:3], Batch = sample_info[, 4]), file = "~/Desktop/Nanopore/Analysis/Differential expression/dCell.csv")
write.csv(x = cbind(sample_info[1:18, 1:3], pca[1:18, 1:3], pca2[, 1:3], Batch = sample_info[1:18, 4]), file = "~/Desktop/Nanopore/Analysis/Differential expression/dCell2.csv")