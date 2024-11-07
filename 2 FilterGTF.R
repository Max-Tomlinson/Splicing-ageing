library(data.table)

cts <- fread(input = "~/Desktop/Nanopore/Bambu/bambu2/counts_transcript.txt")
gtf <- fread(input = "~/Desktop/Nanopore/Bambu/bambu2/extended_annotations.gtf")
gtf$transcript_id <- gsub(pattern = '.*transcript_id "([^"]+)";.*', replacement = '\\1', x = gtf$V9)

expressed_ids <- cts$TXNAME[rowSums(x = cts[, -c(1,2), with = FALSE]) > 0]
expressed_gtf <- gtf[gtf$transcript_id %in% expressed_ids, ]

fwrite(x = expressed_gtf[, 1:9], file = "~/Desktop/Nanopore/Bambu/annotations_bambu.gtf", quote = FALSE, sep = "\t", col.names = FALSE)

novel_expressed_gtf <- expressed_gtf[grepl(pattern = "Bambu", x = transcript_id)]

fwrite(x = novel_expressed_gtf[, 1:9], file = "~/Desktop/Nanopore/Bambu/annotations_bambu_novel.gtf", quote = FALSE, sep = "\t", col.names = FALSE)
