library(dplyr)
library(EnhancedVolcano)
library(gprofiler2)
library(ggrepel)
library(ggsci)

theme_new <- theme_classic(base_size = 24) +
  theme(axis.title.x = element_blank(), 
        legend.position = "bottom", 
        legend.title = element_blank(), 
        text = element_text(family = "Helvetica"))
theme_set(new = theme_new)

results <- NULL

process_file <- function(file) {
  event <- read.table(file = file, header = TRUE)
  title <- tools::toTitleCase(text = gsub(pattern = ".*/bambu_(.*)\\.dpsi", replacement = "\\1", x = file))
  
  dPSI <- colnames(x = event)[1]
  pval <- colnames(x = event)[2]
  
  event$gene_id <- sub(pattern = "\\..*", replacement = "", x = sapply(X = strsplit(x = rownames(x = event), split = ";"), FUN = `[`, 1))
  event <- event[!is.na(x = event[, 1]), ]
  
  gene_ids <- gconvert(query = event$gene_id, organism = "hsapiens", mthreshold = 1, filter_na = TRUE)
  
  event$names <- ifelse(test = !is.na(x = match(x = event$gene_id, table = gene_ids$input)), yes = gene_ids$name[match(x = event$gene_id, table = gene_ids$input)], no = event$gene_id)
  event <- event[order(event[, pval]), ]
  
  sig <- event[event[, pval] <= 0.05, ]
  
  assign(x = "results", value = rbind.data.frame(results, data.frame(
    Event     = title,
    dPSI      = nrow(x = sig),
    dPSI_up   = sum(sig[, dPSI] > 0),
    dPSI_down = sum(sig[, dPSI] < 0)
  )), envir = .GlobalEnv)
  
  if (!startsWith(x = title, prefix = "Isoform")) {
    
    lab_italics <- paste0("italic('", make.unique(names = event$names), "')")
    
    keyvals <- ifelse(test = event[, dPSI] < 0 & event[, pval] < 0.05, yes = "#475da0", no = ifelse(test = event[, dPSI] > 0 & event[, pval] < 0.05, yes = "#ba1f21", no = "grey"))
    
    names(x = keyvals)[keyvals == "#ba1f21"] <- 'Spliced-up'
    names(x = keyvals)[keyvals == "grey"] <- 'NS'
    names(x = keyvals)[keyvals == "#475da0"] <- 'Spliced-down'
    
    plot <- EnhancedVolcano(toptable = event,
                            lab = lab_italics,
                            x = dPSI,
                            y = pval,
                            xlim = range(event[, dPSI], na.rm = TRUE),
                            ylim = c(0, max(-log10(x = event[, pval]), na.rm = TRUE)),
                            xlab = "Difference PSI",
                            selectLab = if (nrow(x = sig) > 20) lab_italics[1:20] else lab_italics[1:nrow(x = sig)],
                            axisLabSize = 24,
                            title = title,
                            titleLabSize = 24,
                            subtitle = NULL,
                            caption = NULL,
                            pCutoff = 0.05,
                            pCutoffCol = pval,
                            FCcutoff = 0,
                            pointSize = c(ifelse(test = abs(x = event[, dPSI]) > 0 & event[, pval] < 0.05, yes = 2, no = 2/8)),
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
    
    ggsave(filename = sub(pattern = ".dpsi", replacement = "_volcano.pdf", x = file), plot = plot, width = 11, height = 8.5, units = "in")
    
  } else {
    
    event <- event[event[, dPSI] > 0, ]
    
    plot <- ggplot(data = event, aes(x = -log10(event[, pval]), y = event[, dPSI])) +
      geom_point(mapping = aes(color = event[, pval] <= 0.05, size = event[, pval] <= 0.05), alpha = 0.75, show.legend = FALSE) +
      geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
      geom_label_repel(mapping = aes(label = ifelse(test = event[, pval] < 0.05, yes = as.character(x = names), no = '')), max.overlaps = 35, force = 10) +
      labs(x = expression(-Log[10]~italic(P)), y = "Differential splicing for transcripts (dPSI)") +
      scale_color_manual(values = c("grey", "#ba1f21")) +
      theme(axis.title.x = element_text())
    
    ggsave(filename = "~/Desktop/Nanopore/Analysis/SUPPA/bambu/bambu_isoform_volcano.pdf", plot = plot, width = 11, height = 8.5, units = "in")
    
  }
}

lapply(X = list.files(path = "~/Desktop/Nanopore/Analysis/SUPPA/bambu", pattern = "dpsi", full.names = TRUE), FUN = process_file)

write.csv(x = results, file = "~/Desktop/Nanopore/Analysis/SUPPA/DSE.csv", row.names = FALSE)

results <- tidyr::pivot_longer(data = results, cols = c(dPSI_up, dPSI_down), names_to = "Direction", values_to = "Count")

ggplot(data = results[!results$Event == "Isoform", ], mapping = aes(x = Event, y = Count)) +
  geom_bar(position = "dodge", stat = "identity", aes(fill = Direction, color = Direction)) +
  labs(x = NULL, y = "Number of differentially spliced events") +
  scale_color_futurama(labels = c("Down", "Up")) +
  scale_fill_futurama(labels = c("Down", "Up"))

ggsave(filename = "~/Desktop/Nanopore/Analysis/SUPPA/DSE.pdf", plot = last_plot(), width = 11, height = 8.5, units = "in")
