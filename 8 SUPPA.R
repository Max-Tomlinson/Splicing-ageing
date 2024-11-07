library(dplyr)
library(ggplot2)
library(ggsci)
library(tidyr)

theme_new <- theme_classic(base_size = 24) +
  theme(axis.title.x = element_blank(), 
        legend.position = "bottom", 
        legend.title = element_blank(), 
        text = element_text(family = "Helvetica"))
theme_set(new = theme_new)

events <- read.table(file = "~/Desktop/Nanopore/Analysis/SUPPA/bambu_all_events.tsv", header = TRUE)
events$Annotation <- ifelse(test = startsWith(x = events$transcript_id, prefix = "ENS"), yes = "Known", no = "Novel")

event_summary <- events %>%
  select(diffsplice_events, gene_id, event_id) %>%
  unique() %>%
  group_by(diffsplice_events) %>%
  summarise(Freq = n()) %>%
  ungroup()

novel_summary <- events %>%
  filter(Annotation == "Novel") %>%
  filter(Type == "inclusion") %>%
  select(diffsplice_events, gene_id, event_id) %>%
  unique() %>%
  group_by(diffsplice_events) %>%
  summarise(Freq = n()) %>%
  ungroup()

data <- merge(event_summary, novel_summary, by = "diffsplice_events", suffixes = c("_known", "_novel"))
data$Known <- data$Freq_known / sum(data$Freq_known)
data$Novel <- data$Freq_novel / sum(data$Freq_novel)
data$Ratio <- data$Freq_novel / data$Freq_known
data <- data %>% gather(Label, Proportion, Known, Novel) %>% mutate(Label = ifelse(test = Label == "Known", yes = "Known", no = "Novel"))

ggplot(data = data, mapping = aes(x = diffsplice_events, y = Proportion)) +
  geom_bar(position = "dodge", stat = "identity", aes(fill = Label, color = Label)) +
  labs(x = NULL, y = "Proportion of alternative splicing events") +
  scale_color_futurama() +
  scale_fill_futurama()

ggsave(filename = "~/Desktop/Nanopore/Analysis/SUPPA/SUPPA_events.pdf", plot = last_plot(), width = 11, height = 8.5)
