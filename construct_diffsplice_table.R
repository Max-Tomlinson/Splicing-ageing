library(tidyverse)
library(magrittr)
library(optparse)

option_list <- list(make_option(c("--suppa_dir"), type = "character", default = NULL, help = "Directory of SUPPA output files without suffix", metavar = "type"))

message(" ## Parsing options")

opt <- optparse::parse_args(object = OptionParser(option_list = option_list))

suppa_dir <- opt$suppa_dir

process_suppa_output <- function(dir, as_event) {
  myas <- read.table(file = paste0(dir, "_" , as_event, "_strict.ioe"), header = TRUE)
  myas$alternative_transcripts <- as.character(x = myas$alternative_transcripts)
  
  s <- strsplit(x = myas$alternative_transcripts, split = ",")
  
  myas <- data.frame(diffsplice_events = as_event,
                     seqname           = rep(x = myas$seqname, sapply(X = s, FUN = length)),
                     gene_id           = rep(x = myas$gene_id, sapply(X = s, FUN = length)),
                     event_id          = rep(x = myas$event_id, sapply(X = s, FUN = length)),
                     total_transcripts = rep(x = myas$total_transcripts, sapply(X = s, FUN = length)),
                     inclusion         = unlist(x = s))
  myas$total_transcripts <- as.character(myas$total_transcripts)
  
  s <- strsplit(x = myas$total_transcripts, split = ",")
  
  myas <- data.frame(diffsplice_events = as_event,
                     seqname           = rep(x = myas$seqname, sapply(X = s, FUN = length)),
                     gene_id           = rep(x = myas$gene_id, sapply(X = s, FUN = length)),
                     event_id          = rep(x = myas$event_id, sapply(X = s, FUN = length)),
                     inclusion         = rep(x = myas$inclusion, sapply(X = s, FUN = length)),
                     exlusion          = unlist(x = s))
  myas$inclusion <- as.character(x = myas$inclusion)
  myas$exlusion  <- as.character(x = myas$exlusion)
  myas <- subset(x = myas, inclusion != exlusion)
  myas <- reshape2::melt(data = myas, measure.vars = c("inclusion", "exlusion"), variable.name = "Type", value.name = "transcript_id")
  myas %<>%
    unique() %>%
    arrange(seqname, gene_id, event_id, transcript_id, Type) %>%
    distinct(event_id, transcript_id, .keep_all = TRUE, fromLast = T) %>%
    dplyr::select(-fromLast)
  }

A3 <- process_suppa_output(dir = suppa_dir, as_event = "A3")
A5 <- process_suppa_output(dir = suppa_dir, as_event = "A5")
AF <- process_suppa_output(dir = suppa_dir, as_event = "AF")
AL <- process_suppa_output(dir = suppa_dir, as_event = "AL")
MX <- process_suppa_output(dir = suppa_dir, as_event = "MX")
RI <- process_suppa_output(dir = suppa_dir, as_event = "RI")
SE <- process_suppa_output(dir = suppa_dir, as_event = "SE")

all_events <- rbind(A3, A5, AF, AL, MX, RI, SE)
write.table(x = all_events, paste0(suppa_dir,"_all_events.tsv"), col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")