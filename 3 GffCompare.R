library(broom)
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggsignif)
library(purrr)
library(readr)
library(tidyr)

theme_new <-  theme_classic(base_size = 24) +
  theme(text = element_text(family = "Helvetica"))
theme_set(new = theme_new)

tmap <- read_tsv(file = "~/Desktop/Nanopore/Bambu/bambu.annotations_bambu.gtf.tmap", show_col_types = FALSE)
tmap <- tmap %>%
  mutate(
    ref_match_len     = if_else(condition = class_code == "u", true = NA_character_, false = ref_match_len),
    length_difference = as.numeric(x = len) - as.numeric(x = ref_match_len),
    qry_gene_type     = if_else(condition = grepl(pattern = "^Bambu", x = qry_gene_id), true = "Novel", false = "Known"),
    qry_type          = if_else(condition = grepl(pattern = "^Bambu", x = qry_id), true = "Novel", false = "Known")
  )

class_code <- tmap %>%
  group_by(class_code, qry_type) %>%
  summarise(
    n              = n(),
    num_exons_mean = mean(num_exons),
    num_exons_sd   = sd(num_exons),
    len_mean       = mean(len),
    len_sd         = sd(len),
    .groups        = 'drop'
  )

length_difference <- tmap %>%
  filter(class_code != "u") %>%
  group_by(class_code, qry_type) %>%
  summarise(
    length_difference_mean = mean(length_difference),
    length_difference_sd   = sd(length_difference),
    .groups                = 'drop'
  )

class_code <- class_code %>%
  left_join(length_difference, by = c("class_code", "qry_type"))

write.csv(x = class_code, file = "~/Desktop/Nanopore/Analysis/GffCompare/tmap.csv", row.names = FALSE)

aov <- aov(formula = log(x = len) ~ class_code, data = tmap)

TukeyHSD <- TukeyHSD(x = aov)
TukeyHSD <- as.data.frame(x = TukeyHSD$`class_code`) %>%
  mutate(sig = ifelse(test = `p adj` < 0.05, yes = "*", no = ""))

write.csv(x = TukeyHSD, file = "~/Desktop/Nanopore/Analysis/GffCompare/TukeyHSD.csv", row.names = TRUE)

ggplot(data = tmap, mapping = aes(x = class_code, y = len, fill = class_code)) +
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.25, fill = "white") +
  labs(x = "Class code", y = "Log transcript length") +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000)) +
  theme(legend.position = "none")

ggsave(filename = "~/Desktop/Nanopore/Analysis/GffCompare/tmap.pdf", plot = last_plot(), width = 11, height = 8.5, units = "in", useDingbats = FALSE)

split <- split(x = tmap, f = tmap$qry_gene_type)

metrics <- c("num_exons", "len", "length_difference")

map1 <- map(.x = metrics, .f = ~ t.test(x = split$Known[[.x]], y = split$Novel[[.x]], alternative = "two.sided")) %>%
  map_dfr(.f = ~ tidy(.x), .id = "metric")

split <- split(x = tmap, f = tmap$qry_type)

map2 <- map(.x = metrics, .f = ~ t.test(x = split$Known[[.x]], y = split$Novel[[.x]], alternative = "two.sided")) %>%
  map_dfr(.f = ~ tidy(.x), .id = "metric")

map <- rbind(map1, map2)
map$metric <- c(paste0(metrics, "_gene"), paste0(metrics, "_transcript"))

write.csv(x = map, file = "~/Desktop/Nanopore/Analysis/GffCompare/map.csv", row.names = FALSE)

output <- read_tsv(file = "~/Desktop/Nanopore/Bambu/output.tab", col_names = c("chr", "start", "end", "strand", "transcript_id", "width", "count"), show_col_types = FALSE)
output <- output %>%
  mutate(
    width = end - start,
    count = stringr::str_count(string = transcript_id, pattern = ",") + 1
  )

chrom_sizes <- read.csv(file = "~/Desktop/Nanopore/References/gene_content.csv")

colnames(x = chrom_sizes) <- c("chr", "size")

output <- left_join(x = output, y = chrom_sizes, by = "chr")

count <- output %>%
  group_by(count) %>%
  summarize(
    n       = sum(count),
    .groups = 'drop'
  )

known <- output %>%
  filter(grepl(pattern = "ENS", x = transcript_id))

write.csv(x = count, file = "~/Desktop/Nanopore/Analysis/GffCompare/count.csv", row.names = FALSE)
write.csv(x = known, file = "~/Desktop/Nanopore/Analysis/GffCompare/known.csv", row.names = FALSE)

chr <- output %>%
  group_by(chr = reorder(chr, match(x = chr, table = c(paste0("chr", 1:22), "chrX")))) %>%
  summarize(n = sum(count),
            n_per_gene = n / unique(x = size) * 100,
            .groups = 'drop') %>%
  mutate(chr = gsub(pattern = "chr", replacement = "", x = chr),
         chr = factor(x = chr, levels = c(1:22, "X")))

ggplot(data = chr, aes(x = chr, y = n)) +
  geom_bar(stat = "identity", colour = "black", fill = "#8da9c6") +
  labs(x = "Chromosome", y = "Number of novel splice junctions")

ggsave(filename = "~/Desktop/Nanopore/Analysis/GffCompare/chr.pdf", plot = last_plot(), width = 11, height = 8.5, units = "in", useDingbats = FALSE)

ggplot(data = chr, aes(x = chr, y = n_per_gene)) +
  geom_bar(stat = "identity", colour = "black", fill = "#8da9c6") +
  labs(x = "Chromosome", y = "Percentage of novel splice junctions per gene")

ggsave(filename = "~/Desktop/Nanopore/Analysis/GffCompare/chr_gene.pdf", plot = last_plot(), width = 11, height = 8.5, units = "in", useDingbats = FALSE)

cpc2 <- read_delim(file = "~/Desktop/Nanopore/Bambu/cpc2_output.txt", delim = "\t", show_col_types = FALSE)
cpc2$label <- factor(x = cpc2$label, levels = unique(x = cpc2$label), labels = c("Non-coding", "Coding"))

ggplot(data = cpc2, aes(x = coding_probability, fill = label)) +
  geom_histogram(bins = 100, colour = "black") +
  labs(x = "Coding probability", y = "Frequency", fill = NULL) +
  scale_fill_futurama() +
  theme(legend.position = "top")

ggsave(filename = "~/Desktop/Nanopore/Analysis/GffCompare/coding_probability.pdf", plot = last_plot(), width = 11, height = 8.5, units = "in", useDingbats = FALSE)

cpc2 <- cpc2 %>%
  gather(key = "type", value = "length", peptide_length, transcript_length) %>%
  mutate(type = factor(type, levels = c("peptide_length", "transcript_length"), labels = c("Protein", "Transcript")))

ggplot(data = cpc2, mapping = aes(x = label, y = length, fill = label)) +
  facet_wrap(~ type, scales = "free_y") +
  geom_boxplot() +
  geom_signif(data = cpc2, stat = "signif", position = "identity", comparisons = list(c("Non-coding","Coding")), map_signif_level = TRUE) +
  labs(x = NULL, y = "Length", fill = NULL) +
  scale_fill_futurama() +
  theme(legend.position = "none")

ggsave(filename = "~/Desktop/Nanopore/Analysis/GffCompare/length.pdf", plot = last_plot(), width = 11, height = 8.5, units = "in", useDingbats = FALSE)
