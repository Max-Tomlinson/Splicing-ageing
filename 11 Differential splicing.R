library(qvalue)

for (file in list.files(path = "~/Desktop/Nanopore/Analysis/glm", pattern = "fitted.txt", full.names = TRUE)) {
  
  results <- data.frame(gene = character(), n = numeric(), diff = numeric(), pval = numeric())
  
  table <- read.table(file = file, header = TRUE, row.names = 1)
  genes <- unique(x = rownames(x = table))
  
  for (gene in genes) {
    
    data <- table[rownames(x = table) == gene, ]
    
    time1 <- as.numeric(x = data[, 1:18])
    time2 <- as.numeric(x = data[, 19:36])
    pairs <- complete.cases(time1, time2)
    
    if (sum(pairs) >= 18) {
      
      diff <- mean(time2 - time1, na.rm = TRUE)
      
      if (diff != 0) {
        
        w <- wilcox.test(x = time1, y = time2, paired = TRUE, exact = FALSE)

        results <- rbind(results, data.frame(genes = gene, Time.2_mean = mean(x = time2, na.rm = TRUE), Time.1_mean = mean(x = time1, na.rm = TRUE), mean_difference = diff, n = sum(pairs), raw_p_values = w$p.value))
      }
    }
  }
  
  results$adjusted_p_values <- qvalue(p = results$raw_p_values)$qvalues
  
  file <- gsub(pattern = ".txt", replacement = "_results.txt", x = file)
  write.table(x = results, file = file, quote = FALSE, sep = "\t", row.names = FALSE)
}